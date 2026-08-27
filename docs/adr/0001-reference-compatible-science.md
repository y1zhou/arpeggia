# Use feature-level scientific references

Arpeggia is reference-compatible: each public scientific calculation names its
own authoritative published method, with an official implementation used to
resolve details the publication leaves unspecified. Deviations require a
recorded rationale; historical behavior and regression tests alone do not make
a deviation intentional. Confirmed scientific corrections may change numerical
or tabular output in a minor release when their effects are documented, while
schema and category changes follow the normal breaking-release policy.

Unresolved scientific concerns are marked `[WARNING]` after maintainer
confirmation when possible. A warning records open debt and never constitutes
acceptance of the behavior.

Reference compatibility is protected by a lean set of curated synthetic and real
fixtures. Goldens record their reference implementation, version or commit,
parameters, and narrow justified tolerance; they are added only when an invariant
or scientific boundary cannot be expressed more directly.
