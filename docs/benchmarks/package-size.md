# Package-size measurements

## v0.9.2 dependency cleanup

These historical measurements compare the v0.9.2 cleanup with
`7e06c82549c1eb94a1e75b2e1d5c2d38f426bc99` (v0.9.1 code), alongside
published v0.9.0/v0.9.1 metadata. They are not v0.10.0 artifact sizes.

Published compressed artifact sizes, in decimal MB:

| Artifact | v0.9.0 | v0.9.1 |
| --- | ---: | ---: |
| CPython 3.13 Linux x86-64 wheel | 6.42 | 21.65 |
| Linux x86-64 GNU CLI archive | 9.37 | 24.65 |

Sources: [v0.9.0 PyPI metadata](https://pypi.org/pypi/arpeggia/0.9.0/json),
[v0.9.1 PyPI metadata](https://pypi.org/pypi/arpeggia/0.9.1/json),
[v0.9.0 GitHub release](https://github.com/y1zhou/arpeggia/releases/tag/v0.9.0),
and [v0.9.1 GitHub release](https://github.com/y1zhou/arpeggia/releases/tag/v0.9.1).

Enabling Polars' `lazy` feature for one NDJSON reader accounted for much of
the avoidable growth. Combined with JSON/Parquet support, it activated query
planning, execution, streaming, and cloud/network support. Removing `lazy`
reduced the normal dependency graph from 311 to 196 package nodes. Parquet,
`kmedoids`, and `sysinfo` remained required for supported features.

The eager replacement preserved schema projection and bounded row-count
validation. An isolated [5B8C benchmark](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/contacts.md#5b8c-contacts-lazy-versus-eager-ndjson)
found essentially unchanged contact-generation timings and unchanged or
faster NDJSON loading on a 2,574-row table.

Same-machine release wheel builds with the same CPython 3.13 interpreter,
Rust toolchain, lockfile, and maturin settings measured:

| Content | Before cleanup | After cleanup |
| --- | ---: | ---: |
| Compressed wheel | 22,376,505 bytes | 9,480,021 bytes |
| Uncompressed extension | 72,859,664 bytes | 34,657,952 bytes |

The wheel was 57.6% smaller. These local wheels target manylinux 2.35; their
absolute sizes are not directly comparable to the published manylinux 2.17
builds. Remaining growth over v0.9.0 is not isolated by this comparison.

At v0.9.2, binary release archives already omitted `docs/` and contained
the binary, README and license. The cleanup extended Cargo exclusions from
benchmark HTML to all `docs/**` and excluded docs from maturin wheels and
source distributions. Archive inspection found no `docs/` entries in the
crate, sdist or wheel, and exactly one extension in the wheel. Source packages
retained structures required by packaged tests and doctests.
