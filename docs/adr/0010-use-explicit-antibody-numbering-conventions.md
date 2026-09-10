# Use explicit antibody numbering conventions

Arpeggia will accept `chothia` as an alias for Martin/enhanced Chothia numbering,
as requested for the antibody-numbering API. This deliberately favors the
structurally corrected convention over compatibility with historical Chothia
outputs. The alias is an Arpeggia API choice: the original conventions remain
distinct in the [scheme authors' numbering service](https://www.bioinf.org.uk/abs/abnum/).

The initial CLI and Python APIs accept named amino-acid strings. A
`NumberedAntibody` represents one variable domain, retaining its input sequence
and domain span. Terminal tags and constant-region tails are allowed;
additional detected variable domains produce an error. Structure-file
integration and constant-region numbering are deferred. Recognizable partial
domains can return results with coverage and diagnostics; numbering itself
does not fill missing sequence.

Germline comparison reports separate V and J reference similarities, coverage
and tied references. It does not reconstruct a unique ancestral antibody or
infer a D segment. References ship offline as a versioned, attributed IMGT
subset covering human and mouse H/K/L and alpaca heavy-chain/VHH references.
Search covers all bundled species by default and accepts an explicit species
restriction. Report matched-reference species rather than presumed input origin.

`cdr_definition="auto"` selects the CDR convention associated with the chosen
numbering scheme. A caller choosing a different CDR definition must explicitly
provide both arguments. This keeps the usual scheme/region pairing convenient
while making a mixed convention intentional. Exact Martin and AHo region
policies still require resolution; the backend's tables are not authoritative
merely because the backend supports those numbering schemes.

An explicit method on `NumberedAntibody` will impute missing residues from the
closest aligned germline reference. Eligible positions, ambiguous reference
matches and the returned result remain under discussion. This method is
separate from ordinary numbering so supplied sequence and inferred sequence
can be distinguished.

Engine selection and the remaining result contract are open. Supporting evidence is in the
[numbering research](https://github.com/y1zhou/arpeggia/blob/master/docs/research/antibody-numbering-schemes-and-tools.md).
