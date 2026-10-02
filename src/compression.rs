//! Streaming compression shared by the FASTA engines and query-file reader.
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use tempfile::NamedTempFile;

use crate::fasta::CancellationToken;

const BUFFER_BYTES: usize = 1 << 20;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Format {
    Plain,
    Gzip,
    Zstd,
}

impl Format {
    pub(crate) fn for_path(path: &Path) -> Self {
        match path
            .extension()
            .and_then(|s| s.to_str())
            .map(str::to_ascii_lowercase)
            .as_deref()
        {
            Some("gz") => Self::Gzip,
            Some("zst" | "zstd") => Self::Zstd,
            _ => Self::Plain,
        }
    }
}

/// Attach the source path to errors that occur after opening the decoder.
struct Input {
    reader: Box<dyn BufRead + Send>,
    path: PathBuf,
}

impl Input {
    fn error(&self, error: io::Error) -> io::Error {
        io::Error::new(
            error.kind(),
            format!("failed to read {}: {error}", self.path.display()),
        )
    }
}

impl Read for Input {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.reader.read(buf).map_err(|error| self.error(error))
    }
}

impl BufRead for Input {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        let path = &self.path;
        self.reader.fill_buf().map_err(|error| {
            io::Error::new(
                error.kind(),
                format!("failed to read {}: {error}", path.display()),
            )
        })
    }

    fn consume(&mut self, amount: usize) {
        self.reader.consume(amount);
    }
}

pub(crate) fn open_input(path: &Path) -> Result<(impl BufRead + Send + use<>, Format)> {
    let mut input = BufReader::with_capacity(
        BUFFER_BYTES,
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    let magic = input
        .fill_buf()
        .with_context(|| format!("failed to read {}", path.display()))?;
    let format = if magic.starts_with(&[0x1f, 0x8b]) {
        Format::Gzip
    } else if magic.starts_with(&[0x28, 0xb5, 0x2f, 0xfd])
        || (magic.len() >= 4 && magic[0] & 0xf0 == 0x50 && magic[1..4] == [0x2a, 0x4d, 0x18])
    {
        Format::Zstd
    } else {
        Format::Plain
    };
    let reader: Box<dyn BufRead + Send> = match format {
        Format::Plain => Box::new(input),
        Format::Gzip => Box::new(BufReader::with_capacity(
            BUFFER_BYTES,
            MultiGzDecoder::new(input),
        )),
        Format::Zstd => Box::new(BufReader::with_capacity(
            BUFFER_BYTES,
            zstd::stream::read::Decoder::with_buffer(input)
                .with_context(|| format!("failed to decode {}", path.display()))?,
        )),
    };
    Ok((
        Input {
            reader,
            path: path.to_owned(),
        },
        format,
    ))
}

enum Encoder {
    Plain(NamedTempFile),
    Gzip(GzEncoder<NamedTempFile>),
    Zstd(zstd::stream::write::Encoder<'static, NamedTempFile>),
}

impl Write for Encoder {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(buf),
            Self::Gzip(writer) => writer.write(buf),
            Self::Zstd(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
            Self::Zstd(writer) => writer.flush(),
        }
    }
}

/// Encoder finalization is explicit: dropping this writer never installs output.
pub(crate) struct AtomicOutput {
    writer: BufWriter<Encoder>,
    destination: PathBuf,
}

impl AtomicOutput {
    pub(crate) fn new(destination: &Path, format: Format) -> Result<Self> {
        let create = || -> Result<Self> {
            let parent = destination
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or_else(|| Path::new("."));
            fs::create_dir_all(parent)?;
            let temporary = NamedTempFile::new_in(parent)?;
            if let Ok(metadata) = fs::metadata(destination) {
                temporary
                    .as_file()
                    .set_permissions(metadata.permissions())?;
            }
            let encoder = match format {
                Format::Plain => Encoder::Plain(temporary),
                Format::Gzip => Encoder::Gzip(GzEncoder::new(temporary, Compression::default())),
                Format::Zstd => Encoder::Zstd(zstd::stream::write::Encoder::new(temporary, 3)?),
            };
            Ok(Self {
                writer: BufWriter::with_capacity(BUFFER_BYTES, encoder),
                destination: destination.to_owned(),
            })
        };
        create().with_context(|| format!("failed to create output {}", destination.display()))
    }

    pub(crate) fn finish(self, cancellation: &CancellationToken) -> Result<()> {
        let destination = self.destination;
        let finish = || -> Result<()> {
            let encoder = self
                .writer
                .into_inner()
                .map_err(|error| error.into_error())?;
            let mut temporary = match encoder {
                Encoder::Plain(temporary) => temporary,
                Encoder::Gzip(encoder) => encoder.finish()?,
                Encoder::Zstd(encoder) => encoder.finish()?,
            };
            temporary.flush()?;
            temporary.as_file().sync_all()?;
            cancellation.check_cancelled()?;
            temporary
                .persist(&destination)
                .map_err(|error| error.error)?;
            Ok(())
        };
        finish().with_context(|| format!("failed to finish output {}", destination.display()))
    }

    fn error(&self, error: io::Error) -> io::Error {
        io::Error::new(
            error.kind(),
            format!("failed to write {}: {error}", self.destination.display()),
        )
    }
}

impl Write for AtomicOutput {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.writer.write(buf).map_err(|error| self.error(error))
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush().map_err(|error| self.error(error))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cancellation_after_writing_never_installs_output() {
        let dir = tempfile::tempdir().unwrap();
        let destination = dir.path().join("output");
        for format in [Format::Plain, Format::Gzip, Format::Zstd] {
            fs::write(&destination, b"existing").unwrap();
            let mut output = AtomicOutput::new(&destination, format).unwrap();
            // Force buffered writes before cancellation, then finalize the encoder.
            output.write_all(&vec![b'A'; BUFFER_BYTES + 1]).unwrap();
            let cancellation = CancellationToken::default();
            cancellation.cancel();
            let error = output.finish(&cancellation).unwrap_err();
            assert!(format!("{error:#}").contains("operation cancelled"));
            assert_eq!(fs::read(&destination).unwrap(), b"existing");
            assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
        }
    }

    #[test]
    fn failed_install_removes_temporary_output_and_reports_destination() {
        let dir = tempfile::tempdir().unwrap();
        let destination = dir.path().join("existing-directory");
        fs::create_dir(&destination).unwrap();
        fs::write(destination.join("sentinel"), b"existing").unwrap();
        for format in [Format::Plain, Format::Gzip, Format::Zstd] {
            let mut output = AtomicOutput::new(&destination, format).unwrap();
            output.write_all(b"new output").unwrap();
            let error = output.finish(&CancellationToken::default()).unwrap_err();
            assert!(
                error
                    .to_string()
                    .contains(&destination.display().to_string())
            );
            assert_eq!(fs::read(destination.join("sentinel")).unwrap(), b"existing");
            assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
        }
    }
}
