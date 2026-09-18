//! Per-connection SQLite scratch routing; never changes process-global settings.
use anyhow::{Context, Result};
use rusqlite::{Connection, OpenFlags, ffi};
use std::{
    ffi::{CString, c_char, c_int, c_void},
    ops::{Deref, DerefMut},
    path::Path,
    ptr,
    sync::{
        Mutex,
        atomic::{AtomicU64, Ordering},
    },
};

static NEXT_VFS: AtomicU64 = AtomicU64::new(0);

struct State {
    parent: *mut ffi::sqlite3_vfs,
    directory: tempfile::TempDir,
    // The underlying OS VFS may retain filename pointers until file close.
    filenames: Mutex<Vec<Box<[u8]>>>,
}

struct ScratchVfs {
    vfs: Box<ffi::sqlite3_vfs>,
    _name: CString,
    #[allow(dead_code)] // Owns the callback state even when diagnostics are disabled.
    state: Box<State>,
}

// Moving the owner never moves its boxed callback state. SQLite's connection
// owns all callback use; filenames are protected for parallel sorter workers.
unsafe impl Send for ScratchVfs {}

impl Drop for ScratchVfs {
    fn drop(&mut self) {
        // ScratchConnection closes SQLite (and sorter workers) before this runs.
        unsafe {
            ffi::sqlite3_vfs_unregister(&mut *self.vfs);
        }
        // TempDir cleans any scratch files still present after file close.
    }
}

unsafe fn state<'a>(vfs: *mut ffi::sqlite3_vfs) -> &'a State {
    unsafe { &*((*vfs).pAppData.cast::<State>()) }
}

unsafe extern "C" fn open(
    vfs: *mut ffi::sqlite3_vfs,
    name: ffi::sqlite3_filename,
    file: *mut ffi::sqlite3_file,
    flags: c_int,
    out: *mut c_int,
) -> c_int {
    // No panic may unwind across SQLite's C ABI.
    std::panic::catch_unwind(|| unsafe {
        let state = state(vfs);
        if !name.is_null() {
            return ((*state.parent).xOpen.unwrap())(state.parent, name, file, flags, out);
        }
        let temporary = match tempfile::Builder::new()
            .prefix("spill-")
            .tempfile_in(state.directory.path())
        {
            Ok(file) => file,
            Err(_) => return ffi::SQLITE_CANTOPEN,
        };
        let path = temporary.path();
        let Some(path) = path.to_str() else {
            return ffi::SQLITE_CANTOPEN;
        };
        // Two terminating NULs meet SQLite's filename/URI contract.
        let mut bytes = path.as_bytes().to_vec();
        bytes.extend_from_slice(&[0, 0]);
        let filename = bytes.into_boxed_slice();
        // Reserve a unique name, then let the OS VFS create it with SQLite's
        // original TEMP/DELETEONCLOSE flags inside our private directory.
        drop(temporary);
        let mut names = match state.filenames.lock() {
            Ok(names) => names,
            Err(_) => return ffi::SQLITE_IOERR,
        };
        names.push(filename);
        let name = names.last().unwrap().as_ptr().cast::<c_char>();
        ((*state.parent).xOpen.unwrap())(state.parent, name, file, flags, out)
    })
    .unwrap_or(ffi::SQLITE_IOERR)
}

// Forward every other callback with the parent's own VFS pointer/pAppData.
macro_rules! forward {
    ($name:ident, $field:ident, ($($arg:ident: $ty:ty),*) -> $ret:ty) => {
        unsafe extern "C" fn $name(vfs: *mut ffi::sqlite3_vfs, $($arg: $ty),*) -> $ret {
            unsafe {
                let parent = state(vfs).parent;
                ((*parent).$field.unwrap())(parent, $($arg),*)
            }
        }
    };
}
type DlSymbol = Option<unsafe extern "C" fn(*mut ffi::sqlite3_vfs, *mut c_void, *const c_char)>;
forward!(delete, xDelete, (name: *const c_char, sync: c_int) -> c_int);
forward!(access, xAccess, (name: *const c_char, flags: c_int, out: *mut c_int) -> c_int);
forward!(full_path, xFullPathname, (name: *const c_char, len: c_int, out: *mut c_char) -> c_int);
forward!(dl_open, xDlOpen, (name: *const c_char) -> *mut c_void);
forward!(dl_error, xDlError, (len: c_int, out: *mut c_char) -> ());
forward!(dl_sym, xDlSym, (handle: *mut c_void, name: *const c_char) -> DlSymbol);
forward!(dl_close, xDlClose, (handle: *mut c_void) -> ());
forward!(random, xRandomness, (len: c_int, out: *mut c_char) -> c_int);
forward!(sleep, xSleep, (micros: c_int) -> c_int);
forward!(time, xCurrentTime, (out: *mut f64) -> c_int);
forward!(error, xGetLastError, (len: c_int, out: *mut c_char) -> c_int);
forward!(time64, xCurrentTimeInt64, (out: *mut i64) -> c_int);
forward!(set_call, xSetSystemCall, (name: *const c_char, call: ffi::sqlite3_syscall_ptr) -> c_int);
forward!(get_call, xGetSystemCall, (name: *const c_char) -> ffi::sqlite3_syscall_ptr);
forward!(next_call, xNextSystemCall, (name: *const c_char) -> *const c_char);

/// Field order is intentional: close the connection before unregistering VFS
/// callbacks and deleting its private scratch directory.
pub(crate) struct ScratchConnection {
    connection: Option<Connection>,
    _vfs: ScratchVfs,
}
impl Deref for ScratchConnection {
    type Target = Connection;
    fn deref(&self) -> &Connection {
        self.connection.as_ref().unwrap()
    }
}
impl DerefMut for ScratchConnection {
    fn deref_mut(&mut self) -> &mut Connection {
        self.connection.as_mut().unwrap()
    }
}
impl ScratchConnection {
    pub(crate) fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let folder = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let folder = std::fs::canonicalize(folder).with_context(|| {
            format!("cannot resolve SQLite scratch folder {}", folder.display())
        })?;
        let directory = tempfile::Builder::new()
            .prefix(".taxutils-sqlite-")
            .tempdir_in(&folder)
            .with_context(|| {
                format!(
                    "cannot create SQLite scratch directory in {}",
                    folder.display()
                )
            })?;
        let name = CString::new(format!(
            "taxutils-scratch-{}",
            NEXT_VFS.fetch_add(1, Ordering::Relaxed)
        ))?;
        unsafe {
            anyhow::ensure!(
                ffi::sqlite3_initialize() == ffi::SQLITE_OK,
                "SQLite initialization failed"
            );
            let parent = ffi::sqlite3_vfs_find(ptr::null());
            anyhow::ensure!(!parent.is_null(), "SQLite has no default VFS");
            let mut state = Box::new(State {
                parent,
                directory,
                filenames: Mutex::new(Vec::new()),
            });
            let mut vfs = Box::new(ptr::read(parent));
            vfs.pNext = ptr::null_mut();
            vfs.zName = name.as_ptr();
            vfs.pAppData = (&mut *state as *mut State).cast();
            vfs.xOpen = Some(open);
            macro_rules! bind {
                ($field:ident,$callback:ident) => {
                    if vfs.$field.is_some() {
                        vfs.$field = Some($callback);
                    }
                };
            }
            bind!(xDelete, delete);
            bind!(xAccess, access);
            bind!(xFullPathname, full_path);
            bind!(xDlOpen, dl_open);
            bind!(xDlError, dl_error);
            bind!(xDlSym, dl_sym);
            bind!(xDlClose, dl_close);
            bind!(xRandomness, random);
            bind!(xSleep, sleep);
            bind!(xCurrentTime, time);
            bind!(xGetLastError, error);
            bind!(xCurrentTimeInt64, time64);
            bind!(xSetSystemCall, set_call);
            bind!(xGetSystemCall, get_call);
            bind!(xNextSystemCall, next_call);
            anyhow::ensure!(
                ffi::sqlite3_vfs_register(&mut *vfs, 0) == ffi::SQLITE_OK,
                "cannot register scratch VFS"
            );
            let owner = ScratchVfs {
                vfs,
                _name: name,
                state,
            };
            let connection = Connection::open_with_flags_and_vfs(
                path,
                OpenFlags::default(),
                owner._name.to_str()?,
            )
            .with_context(|| format!("cannot open {} with local SQLite scratch", path.display()))?;
            Ok(Self {
                connection: Some(connection),
                _vfs: owner,
            })
        }
    }

    pub(crate) fn close(mut self) -> rusqlite::Result<()> {
        match self.connection.take().unwrap().close() {
            Ok(()) => Ok(()),
            Err((connection, error)) => {
                self.connection = Some(connection);
                Err(error)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spill(connection: &ScratchConnection) {
        connection
            .execute_batch(
                "PRAGMA temp_store=FILE; PRAGMA cache_size=-64;
             PRAGMA temp.cache_size=-64; PRAGMA threads=2;
             CREATE TEMP TABLE stage(k INTEGER, payload BLOB);
             WITH RECURSIVE seq(x) AS (VALUES(1) UNION ALL SELECT x+1 FROM seq WHERE x<60000)
             INSERT INTO stage SELECT x,randomblob(256) FROM seq;",
            )
            .unwrap();
        let staged = connection._vfs.state.filenames.lock().unwrap().len();
        assert!(staged > 0, "TEMP table did not spill");
        connection
            .execute_batch(
                "CREATE TABLE data AS SELECT * FROM stage;
             CREATE INDEX idx_payload ON data(payload);",
            )
            .unwrap();
        assert!(
            connection._vfs.state.filenames.lock().unwrap().len() > staged,
            "index sorter did not create scratch files"
        );
        let root = connection._vfs.state.directory.path();
        for name in connection._vfs.state.filenames.lock().unwrap().iter() {
            let path = std::ffi::CStr::from_bytes_until_nul(name)
                .unwrap()
                .to_str()
                .unwrap();
            assert!(Path::new(path).starts_with(root));
        }
        #[cfg(target_os = "linux")]
        {
            // Verify OS-level descriptors, including SQLite's unlinked TEMP file.
            let targets: Vec<_> = std::fs::read_dir("/proc/self/fd")
                .unwrap()
                .filter_map(|f| std::fs::read_link(f.ok()?.path()).ok())
                .collect();
            assert!(
                targets.iter().any(|p| p.starts_with(root)),
                "no actual scratch descriptor in selected directory"
            );
        }
        assert_eq!(
            connection
                .query_row("PRAGMA integrity_check", [], |r| r.get::<_, String>(0))
                .unwrap(),
            "ok"
        );
    }

    #[test]
    fn real_temp_table_and_index_spills_stay_local_and_cleanup() {
        let root = tempfile::tempdir().unwrap();
        let connection = ScratchConnection::open(root.path().join("db.sqlite")).unwrap();
        let scratch = connection._vfs.state.directory.path().to_path_buf();
        spill(&connection);
        connection.close().unwrap();
        assert!(!scratch.exists());
        assert!(root.path().join("db.sqlite").exists());
    }

    #[test]
    fn concurrent_connections_keep_independent_directories() {
        let roots: Vec<_> = (0..2).map(|_| tempfile::tempdir().unwrap()).collect();
        let barrier = std::sync::Arc::new(std::sync::Barrier::new(2));
        let threads: Vec<_> = roots
            .iter()
            .map(|root| {
                let folder = root.path().to_path_buf();
                let barrier = barrier.clone();
                std::thread::spawn(move || {
                    let connection = ScratchConnection::open(folder.join("db.sqlite")).unwrap();
                    barrier.wait();
                    spill(&connection);
                })
            })
            .collect();
        for thread in threads {
            thread.join().unwrap();
        }
        for root in roots {
            assert!(std::fs::read_dir(root.path()).unwrap().all(|e| {
                !e.unwrap()
                    .file_name()
                    .to_string_lossy()
                    .starts_with(".taxutils-sqlite-")
            }));
        }
    }

    #[test]
    fn unavailable_scratch_errors_without_falling_back() {
        let root = tempfile::tempdir().unwrap();
        let connection = ScratchConnection::open(root.path().join("db.sqlite")).unwrap();
        std::fs::remove_dir(connection._vfs.state.directory.path()).unwrap();
        let error = connection
            .execute_batch(
                "PRAGMA temp_store=FILE; PRAGMA temp.cache_size=-16;
             CREATE TEMP TABLE stage(payload);
             WITH RECURSIVE seq(x) AS (VALUES(1) UNION ALL SELECT x+1 FROM seq WHERE x<10000)
             INSERT INTO stage SELECT randomblob(1024) FROM seq;",
            )
            .unwrap_err();
        assert_eq!(
            error.sqlite_error_code(),
            Some(rusqlite::ErrorCode::CannotOpen)
        );
        assert!(connection._vfs.state.filenames.lock().unwrap().is_empty());
    }
}
