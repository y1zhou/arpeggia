# Antibody numbering engine comparison

On 10 September 2026, AntPack, Immunum and RIOT were compared on the same
26,365 protein inputs. Immunum had the highest agreement with the supplied
numbering labels in every tested scheme. RIOT rejected all negative controls.
Neither result establishes an independent numbering-accuracy ranking.

## Versions and settings

| Engine | Tested version | Configuration |
| --- | --- | --- |
| [AntPack](https://pypi.org/project/antpack/0.3.8.6.3/) | 0.3.8.6.3, CPython 3.13.13 | `SingleChainAnnotator(chains=['H', 'K', 'L'], scheme=...)`; `analyze_seq()` |
| [Immunum](https://github.com/ENPICOM/immunum/tree/45bb70d34802cc592ebd86e685cc9f551885a2d6) | 1.3.1, native Rust release build | Default features disabled; H/K/L; default confidence threshold 0.5 |
| [RIOT](https://pypi.org/project/riot-na/5.2.0/) | 5.2.0, CPython 3.12.13 | Protein API; bundled human/mouse/alpaca references; `extend_alignment=True`, `return_all_domains=False` |

AntPack's comparator is the GPL release requiring no license key; results do
not describe its newer releases. RIOT's inspected source declares 5.3.0, but
5.2.0 was the latest published package at retrieval. Its numbering implementation
and 68 database files match that source; its alignment backend differs. AHo is
unsupported by RIOT. Dependencies and benchmark adapters were isolated from
Arpeggia's environment.

## Inputs and comparison rules

Fixtures are pinned to AntPack commit
[`a136d1f`](https://github.com/jlparkI/AntPack/tree/a136d1f066400edc110bb6b30f82b194f7cb4c3a/tests/test_data).
All source rows are retained: 24,083 distinct sequences across 26,365 rows.
The normalized JSONL input SHA-256 is
`f0cc87c76f53801355b83d16e16c25bfd9d91378eb22b5d74e2fcbaba1689a2a`.

| Fixture group | Rows | Source within the pinned fixture directory |
| --- | ---: | --- |
| Numbering labels | 1,484 | `test_data.csv.gz`, `sequence` column |
| Named therapeutic chains | 1,190 | `addtnl_test_data.fasta.gz` |
| V/J examples | 892 | `vj_gene_testing.csv.gz`, `sequence` column |
| COVID antibody chains | 20,680 | `covid_data.fasta.gz` |
| Scoring examples | 1,148 | `imgt_comp_scoring.csv.gz`, `sequences` column |
| TCR controls | 525 | `tcr_test_data.csv.gz`, `sequence` column |
| Other-protein controls | 356 | `non_antibody_test_data/converted_seqs.txt.gz`; remove alignment gaps |
| Terminal truncations | 90 | Derived from numbering fixtures as described below |

Truncation parents are the first ten H, K and L sequences each with fixture
IMGT end labels 1 and 127/128; AntPack supplies chain classification. Remove
five residues from the beginning, end, or both. DNA translation and TCR-specific
germline/CDR fixtures are outside this protein-antibody comparison. Prepared
clustering sequences are not added again. Groups overlap substantially, so
their row counts must not be treated as independent biological samples.

Each result becomes one position or null per original input residue. Inclusive
Immunum spans become half-open spans; RIOT input offsets and mapped amino acids
are checked against the original sequence. Decimal and letter insertions are
equivalent (`111.1 = 111A`, `111.27 = 111AA`). Compare in input order, preserving
IMGT's reversed insertion order around position 112. A returned mapping means
at least one numbered residue. Advisory messages are retained separately;
discarding warned outputs would change the sample being compared.

The [upstream numbering test](https://github.com/jlparkI/AntPack/blob/a136d1f066400edc110bb6b30f82b194f7cb4c3a/tests/numbering_tools/test_single_chain_annotator.py)
describes PDB sequences numbered by another tool without identifying its
version. These are comparison labels, not independently adjudicated truth.
Numbering row 671 contains 80 residues but 81 IMGT and Kabat labels. Exclude
only those two malformed label arrays from agreement denominators; retain the
sequence in all engine runs.

## Agreement with supplied numbering labels

Exact agreement requires every residue label, including unnumbered residues,
to match. Parentheses give the percentage of valid fixture rows.

| Scheme | AntPack | Immunum | RIOT |
| --- | ---: | ---: | ---: |
| IMGT | 1,444/1,483 (97.37%) | 1,458/1,483 (98.31%) | 1,448/1,483 (97.64%) |
| Martin | 1,448/1,484 (97.57%) | 1,459/1,484 (98.32%) | 1,391/1,484 (93.73%) |
| Kabat | 1,443/1,483 (97.30%) | 1,460/1,483 (98.45%) | 1,392/1,483 (93.86%) |
| AHo | 1,462/1,484 (98.52%) | 1,478/1,484 (99.60%) | Unsupported |

| Scheme | Compared residues | AntPack label agreement | Immunum label agreement | RIOT label agreement |
| --- | ---: | ---: | ---: | ---: |
| IMGT | 170,015 | 99.871% | 99.954% | 99.920% |
| Martin | 170,095 | 99.900% | 99.905% | 99.800% |
| Kabat | 170,015 | 99.909% | 99.951% | 99.797% |
| AHo | 170,095 | 99.958% | 99.983% | Unsupported |

## Recognition, coverage and failures

| Returned mappings | AntPack | Immunum | RIOT |
| --- | ---: | ---: | ---: |
| Antibody/fragment inputs, IMGT | 25,484/25,484 | 25,484/25,484 | 25,484/25,484 |
| TCR controls, IMGT | 525/525 | 0/525 | 0/525 |
| Other-protein controls, IMGT | 356/356 | 42/356 | 0/356 |

Negative-control counts are unchanged across supported schemes. All 90
truncations preserve each engine's parent numbering and chain under every
supported scheme. The three engines also agree with each other on all 90 IMGT
truncations. This establishes consistency for modest terminal losses, not
correct reconstruction of missing sequence.

Immunum returns nine inconsistent Martin results and nine inconsistent Kabat
results on the COVID set: its reported domain contains one more residue than
its position array. The adapter rejects these mappings, leaving 25,475 positive
results per scheme. A direct native-library probe confirms an upstream defect.
For COVID row 17116 (`6613_L`), the 107-residue input ending `FGPGTKVDIKR`
receives only 106 labels while its inclusive span remains `0..=106`. The light-chain
conversion tables stop at source IMGT position 127 and omit aligned position
128. All nine cases share this pattern. AntPack/RIOT leave that `R` unnumbered
in Martin/Kabat; AntPack and Immunum include it as AHo 149. Removing the terminal
`R` makes all four Immunum schemes consistent. A correction can preserve it as
an excluded tail after verifying that only terminal positions were dropped;
blindly shortening spans to the label count could conceal internal omissions.
The relevant upstream code is the
[Martin light-chain rules](https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/numbering/martin.rs#L150)
and [result construction](https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/annotator.rs#L155).

All 42 Immunum non-antibody hits number only 1–6 residues near the input's end;
18 report confidence 1.0. For example, other-protein row 4 maps only the final
`D` of a 172-residue input to light-chain position 1. Their confidence range
0.5154–1.0 overlaps positive inputs (0.5512–0.9683). The
[normalized score](https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/alignment.rs#L260-L299)
measures the selected alignment without requiring meaningful domain coverage.
Raising its cutoff cannot separate these cases; input length validation also
does not constrain the returned domain's length.

AntPack returns mappings for every negative input in this H/K/L configuration;
381 of the 525 TCR results have no IMGT advisory message. All negative controls
have template identity below 0.85, but imposing that cutoff would also reject
10 of the 90 valid truncations. Neither an empty warning nor a cutoff tuned on
these fixtures establishes a reliable antibody classifier.

RIOT reports no exceptions, but some returned results carry validation flags.
Its 634 locus-name flags per scheme arise from alpaca reference names lacking
the literal `IGH`; 845 containment flags compare extended numbering with the
unextended alignment. These flags do not establish incorrect residue mapping.
Insertion-location flags occur in 99 IMGT, 155 Martin and 424 Kabat results and
merit inspection alongside positional disagreements.

## Disagreements with downstream consequences

These IMGT examples use one-based fixture row numbers. No engine is treated as
truth, and majority agreement does not resolve the discrepancy.

| Therapeutic row | Observed disagreement | Consequence |
| --- | --- | --- |
| 200, frexalimab heavy | Input residue 34 is 38 in Immunum, 39 in AntPack/RIOT; all place conserved W41 at input residue 36 | A residue changes CDR1/FR2 membership despite identical domain coverage |
| 749, elipovimab light | Input residues 61–62 are 82A/83 in AntPack, 83/83A in Immunum, 80/83 in RIOT | Framework insertion placement differs in all three; RIOT flags an earlier insertion |
| 936, nimotuzumab light | AntPack/RIOT number terminal T as 127; Immunum leaves it outside the domain | FR4 coverage and subsequent imputation can differ |

For the 1,190 therapeutic inputs, exact IMGT maps agree in 1,179 AntPack/Immunum,
1,177 AntPack/RIOT and 1,173 Immunum/RIOT comparisons. Coverage disagreements
occur in 7, 3 and 8 cases respectively; equal labels over shared residues do not
imply identical domain boundaries.

An exact-sequence join to the fixture clustering metadata identifies 581 COVID
rows associated only with alpaca heavy-V references, excluding ambiguous
species associations. Exact IMGT maps agree in 572 Immunum/AntPack, 561
Immunum/RIOT and 559 AntPack/RIOT comparisons; all return heavy-chain mappings.
These are reference-species annotations, not independently verified input
organisms or structural numbering labels.

## Recognition follow-up and partial loops

On 11 September, a native probe repeated all 26,365 inputs through Immunum's
public raw-alignment API. It reproduced every confidence decision and the
accepted chains/spans. All 25,484 positive inputs matched 80–121 distinct
profile positions; the 42 false hits matched only 1–6. Requiring 30 matched
positions alongside confidence 0.5 rejected all 881 negative controls without
losing an original positive input. Any threshold from 7 through 80 separated
this panel, so it does not establish an optimal threshold. Insertions do not
count as matched profile positions.

For each of the earlier 30 H/K/L parents, retain the first or last 20, 30, 40,
60 or 80 residues, producing 300 additional fragments. Compare safely converted
IMGT labels with each parent's cropped numbering, including unnumbered positions.
These are consistency references, not independent accuracy labels. The table
uses the tentative 30-position gate plus confidence 0.5; each group has 30 inputs.
The 60 diagnostic 20-residue inputs fall below the high-level engine's minimum
input length and fail the combined acceptance check.

| Retained length | First residues: pass gates | Preserve parent labels | Last residues: pass gates | Preserve parent labels |
| --- | ---: | ---: | ---: | ---: |
| 30 | 30/30 | 0/30 | 21/30 | 21/21 |
| 40 | 30/30 | 30/30 | 29/30 | 29/29 |
| 60 | 29/30 | 25/29 | 30/30 | 25/30 |
| 80 | 30/30 | 30/30 | 30/30 | 11/30 |

All fragments retain their parent chain classification. Among the 229 passing
fragments, only 171 preserve final numbering. Of these passing fragments, 28
have unchanged raw alignment states but changed final labels: the length-based
converter can treat a terminally cut CDR as a complete shorter loop. For
numbering parent 10's first 30 residues, surviving IMGT 30/31 become 37/38.

Using IMGT CDR ranges 27–38, 56–65 and 105–117, only 7/64 passing fragments
with a raw endpoint inside a CDR preserve parent labels. Framework endpoints
fare better, at 164/165, but are insufficient: parent 12's first 60 residues
have raw endpoints 1–66 while their last five residues change from parent
61–65 to 62–66. The inferred endpoint has moved just into FR3.

A coverage gate rejects trivial matches but does not resolve missing loop
context. Partial-domain support needs an explicit scope and diagnostic policy;
neither confidence nor endpoint classification guarantees stable numbering.
The earlier modest FR1/FR4 truncations remain consistent in all 90 cases.
Requiring raw coverage from FR1 through FR4 (`cons_start <= 26`,
`cons_end >= 118`) would also reject ten original antibody-panel rows. This is
a possible scope restriction, not an accepted policy or proof of correctness.

## Implications for Arpeggia

Immunum's Rust core was selected after this comparison: four requested schemes
and strongest fixture-label agreement. Integration still requires resolving its
conversion/span defect, requiring meaningful aligned-domain coverage,
and addressing the long-insertion and multiple-domain limitations documented
in the [numbering research](https://github.com/y1zhou/arpeggia/blob/master/docs/research/antibody-numbering-schemes-and-tools.md).
RIOT is a useful comparator with stronger negative-control rejection here;
this panel does not show better positional accuracy than Immunum.

The panel has no independent species-specific structural truth. Immunum's
positive domain spans are 80–139 residues; these inputs do not qualify
extreme loop lengths, severe partial domains, multidomain rejection,
germline matching or imputation. Runs used different interfaces and concurrent
processes without repeated timing, so no comparative speed claim is made.
Arpeggia owns the bindings and rendering under
[ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md).
