//! Central thread-count policy.
//!
//! Every parallel stage in this crate takes its worker count from here, so the
//! number of threads a run uses is a single value the caller chooses rather
//! than a per-site heuristic. `None` means "all logical CPUs".
//!
//! Pools are built per call and used through [`install`]; the crate never
//! installs a global Rayon pool, so a thread count applies to the operation it
//! was passed to and nothing else in the host process.

use anyhow::{Result, bail};

/// Turn a caller's request into a concrete worker count.
///
/// `None` resolves to the machine's logical CPU count, falling back to 1 when
/// the platform cannot report it.
pub fn resolve(threads: Option<usize>) -> Result<usize> {
    match threads {
        Some(0) => bail!("threads must be at least 1"),
        Some(value) => Ok(value),
        None => Ok(std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1)),
    }
}

/// Run `operation` on a private Rayon pool of `threads` workers.
///
/// Using a scoped pool rather than `build_global` keeps the choice local to
/// this call and makes the result independent of `RAYON_NUM_THREADS`.
pub fn install<T, F>(threads: usize, operation: F) -> Result<T>
where
    F: FnOnce() -> T + Send,
    T: Send,
{
    Ok(rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?
        .install(operation))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_honors_explicit_counts_and_rejects_zero() {
        assert_eq!(resolve(Some(3)).unwrap(), 3);
        assert!(resolve(Some(0)).is_err());
        assert!(resolve(None).unwrap() >= 1);
    }

    #[test]
    fn install_runs_on_a_pool_of_the_requested_width() {
        let width = install(3, rayon::current_num_threads).unwrap();
        assert_eq!(width, 3);
    }
}
