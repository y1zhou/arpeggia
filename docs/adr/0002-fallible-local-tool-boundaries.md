# Use fallible APIs for a local scientific tool

Arpeggia is primarily a local scientific CLI and library, not an intrinsically
multi-tenant execution sandbox. Public Rust APIs return typed errors for invalid
input instead of panicking, and the CLI and Python bindings translate those
errors into their native error forms. The core validates finite, physically
meaningful parameter domains but does not impose arbitrary workload ceilings;
services embedding Arpeggia own quotas, timeouts, input-size limits, and process
resource controls.

Changing existing Rust signatures is accepted now rather than preserving
panic-based compatibility, because predictable scientific failure modes take
priority over the existing error contract.
