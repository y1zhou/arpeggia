# Use fallible APIs for a local scientific tool

Arpeggia is primarily a local scientific CLI and library, not an intrinsically
multi-tenant execution sandbox. Public Rust APIs return typed errors for invalid
input instead of panicking, and the CLI and Python bindings translate those
errors into their native error forms. The core validates finite, physically
meaningful parameter domains but does not impose arbitrary workload ceilings;
services embedding Arpeggia own quotas, timeouts, input-size limits, and process
resource controls. Algorithm-specific supported-input bounds, such as the
[antibody domain limits](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md#inputs-and-result-objects),
are separate from service quotas.

A calculation that cannot produce a complete scientifically meaningful scalar
is a typed `Calculation` failure. Python raises `RuntimeError` and the CLI exits
nonzero; neither boundary substitutes `None` or a partial score. Python
`ValueError` remains reserved for invalid public arguments.

Predictable scientific failures take precedence over compatibility with
panic-based Rust signatures.
