# Build and test

Install Rust and [uv](https://docs.astral.sh/uv/getting-started/installation/),
then clone the repository. Checked-in lock files define dependency versions.

```bash
git clone https://github.com/y1zhou/arpeggia.git
cd arpeggia
cargo build --locked
cargo test --locked
uv sync --frozen --all-extras
uv run --no-sync maturin develop --uv --features python --locked
uv run --no-sync pytest python/tests
uv run --no-sync ty check python
```

Rerun `maturin develop` after changing Rust or switching branches with native API
changes. It rebuilds the extension used by `uv run --no-sync python`; an editable
Python install alone does not refresh native code.

For optimized builds, add `--release` to `maturin develop`, or build artifacts:

```bash
cargo build --release --locked
uv run --no-sync maturin build --release --features python --locked
```

Install the CLI from this checkout with `cargo install --path . --locked`.
See the [CI workflows](https://github.com/y1zhou/arpeggia/tree/master/.github/workflows)
for release targets and the [README](https://github.com/y1zhou/arpeggia/blob/master/README.md#usage)
for usage.
