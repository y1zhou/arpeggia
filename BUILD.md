# Build Instructions

This document provides detailed instructions for building arpeggia as both a Rust library/CLI and a Python package.

## Prerequisites

- Rust toolchain (1.91+)
- Python 3.10+ (for Python bindings)
- Python development headers (for Python bindings)

### Installing Python Development Headers

**Ubuntu/Debian:**
```bash
sudo apt-get install python3-dev
```

**macOS:**
```bash
brew install python
```

**Windows:**
Python development headers are included with the standard Python installer.

## Building the Rust CLI/Library

### Standard Build (CLI only)

```bash
cargo build --release
```

This builds the `arpeggia` CLI binary without Python support. The binary will be located at `target/release/arpeggia`.

### Install CLI Globally

```bash
cargo install --path .
```

This installs the CLI to `~/.cargo/bin/arpeggia`.

### Running Tests

```bash
cargo test
```

## Building the Python Package

### Option 1: Development Build with Maturin

This is the recommended approach for development:

```bash
# Install maturin
pip install maturin

# Build and install in development mode
maturin develop --release --features python
```

This builds the package and installs it in your current Python environment in editable mode.

### Option 2: Build Wheel

To create a distributable wheel file:

```bash
maturin build --release --features python
```

The wheel will be created in `target/wheels/` and can be installed with:

```bash
pip install target/wheels/arpeggia-*.whl
```

### Option 3: Build and Install in One Step

```bash
maturin develop --release --features python
```

## Testing the Python Package

After building with maturin, run the test script:

```bash
python python/test_arpeggia.py
```

Or test manually in Python:

```python
import arpeggia
import polars as pl

# Test with a PDB file
contacts = arpeggia.contacts("test-data/1ubq.pdb")
print(f"Found {len(contacts)} contacts")
print(contacts.head())
```

## Common Build Issues

### Issue: `pyo3` or `pyo3-polars` version conflicts

**Solution:** The `Cargo.toml` is configured with compatible versions. If you encounter issues:
- Ensure you're using polars 0.48.0
- Ensure pyo3 0.24 and pyo3-polars 0.21 are being used
- Run `cargo clean` and rebuild

### Issue: Python headers not found

**Error:** `fatal error: Python.h: No such file or directory`

**Solution:** Install Python development headers (see Prerequisites above).

### Issue: `cargo build` with `--features python` fails

**Solution:** Make sure Python 3 is in your PATH:
```bash
which python3  # Should show Python location
python3 --version  # Should show Python 3.10+
```

## Version Compatibility

| Component | Version | Notes |
|-----------|---------|-------|
| Rust | 1.91+ | Minimum required version |
| Python | 3.10+ | For Python bindings |
| Polars | 0.52.0 | Latest stable version |
| PyO3 | 0.26 | Latest stable version |
| pyo3-polars | 0.25 | Matches polars 0.52 |

## Continuous Integration

For CI/CD pipelines:

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Build CLI
cargo build --release --locked

# Build Python wheel (requires maturin)
pip install maturin
maturin build --release --features python

# Run tests
cargo test
```

## Cross-Compilation

For cross-platform wheels, use `maturin` with appropriate targets:

```bash
# Build for multiple platforms
maturin build --release --features python --target x86_64-unknown-linux-gnu
maturin build --release --features python --target aarch64-apple-darwin
```

## Troubleshooting

If you encounter build issues:

1. **Clean build artifacts:**
   ```bash
   cargo clean
   rm -rf target/
   ```

2. **Update dependencies:**
   ```bash
   cargo update
   ```

3. **Verify Python environment:**
   ```bash
   python3 -c "import sysconfig; print(sysconfig.get_config_var('INCLUDEPY'))"
   ```

4. **Check cargo tree for conflicts:**
   ```bash
   cargo tree | grep polars
   cargo tree | grep pyo3
   ```

## Publishing

### Publishing to crates.io (Rust)

```bash
cargo publish
```

### Publishing to PyPI (Python)

```bash
# Build wheels for multiple platforms
maturin build --release --features python

# Upload to PyPI
maturin publish --features python
```
