use crate::ArpeggiaError;
use polars::prelude::*;

pub(crate) fn string_values<'a>(values: impl Iterator<Item = &'a str>) -> StringChunked {
    StringChunked::from_iter_values(PlSmallStr::EMPTY, values)
}

pub(crate) fn polars_calculation_error(error: PolarsError) -> ArpeggiaError {
    ArpeggiaError::Calculation(error.to_string())
}

/// Execute a parallel operation with the configured thread limit.
/// Uses rayon's thread pool with the specified number of threads.
pub fn run_with_threads<F, T>(num_threads: isize, f: F) -> T
where
    F: FnOnce() -> T + Send,
    T: Send,
{
    let n_threads = num_threads.max(0) as usize;
    if n_threads == 0 || n_threads >= rayon::current_num_threads() {
        // Use global pool directly when auto or enough threads
        f()
    } else {
        // Create a scoped pool with limited threads
        match rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build()
        {
            Ok(pool) => pool.install(f),
            Err(_) => f(),
        }
    }
}
