# Build and test

The checked-in lock files define release dependency state.

```bash
cargo build --locked
cargo test --locked
uv sync --frozen --all-extras
uv run maturin develop --uv --features python --locked
uv run --no-sync pytest python/tests
```

Release builds use the same locked inputs:

```bash
cargo build --release --locked
uv run maturin build --release --features python --locked
```

See `README.md` for command usage and `.github/workflows/` for the supported
release targets.
