# Antibody numbering engine comparison

On 10 September 2026, AntPack, Immunum and RIOT were compared on the same
26,365 protein inputs. Immunum had the highest agreement with the supplied
numbering labels in every tested scheme. RIOT rejected all negative controls.
Neither result establishes an independent numbering-accuracy ranking.

## Chothia adapter qualification

On 16 September 2026, the release Python wheel numbered the same pinned 26,365
inputs with distinct `scheme="chothia"`, automatic CDRs and germline matching.
It accepted 25,474 inputs and rejected all 881 negative controls plus ten
fragments lacking IMGT 23–118 anchor coverage. Every accepted result retained
contiguous input correspondence, unique ordered labels and valid region names;
all had a V match and five lacked qualifying J evidence.

Regression checks distinguish heavy H82 insertions from Martin H72 insertions
and the schemes' short light-CDR1 deletion ordering. They also cover mixed CDR
definitions, the light-chain terminal-span correction, Chothia's 26-insertion
limit, terminal imputation and rejection of mixed-scheme antibody alignments.
These qualify adapter behavior; Chothia was not part of the original
three-engine positional-agreement comparison below.

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
the initial scope restriction, later tightened to IMGT 23–118 in the
[terminal-framework follow-up](#terminal-framework-regression-follow-up). Neither
boundary proves numbering correctness.

## Implications for Arpeggia

Immunum's Rust core was selected after this comparison: four requested schemes
and strongest fixture-label agreement. The adapter qualification below covers
its conversion/span correction,
aligned-domain coverage, insertion limits and multiple-domain guards.
The [numbering research](https://github.com/y1zhou/arpeggia/blob/master/docs/research/antibody-numbering-schemes-and-tools.md)
records engine evidence; [ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md)
records accepted scope. RIOT is a useful comparator with stronger
negative-control rejection here;
this panel does not show better positional accuracy than Immunum.

The panel has no independent species-specific structural truth. Immunum's
positive domain spans are 80–139 residues; these inputs do not qualify
extreme loop lengths, severe partial domains, multidomain rejection,
germline matching or imputation. Runs used different interfaces and concurrent
processes without repeated timing, so no comparative speed claim is made.
Arpeggia owns the bindings and rendering under
[ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md).

## Arpeggia adapter qualification

On 11 September 2026, Arpeggia's release Python wheel repeated all 26,365 inputs
in all four schemes: 105,460 calls, retaining errors, diagnostics and one position
or null per input residue. This used the same pinned inputs and normalization
above, default CDR conventions, and the then-bundled human, mouse and alpaca
reference species. Numbering and
matching were validated at implementation milestone `5b798a3`; the subsequent
`a9f1337` change adds an imputation coverage diagnostic without changing numbering.

Every scheme accepts 25,474 antibody inputs and rejects all 881 negative
controls. The other ten rejections are the agreed severe-partial restriction:
numbering rows 398, 599, 613, 671 and 1077; therapeutic row 1107; and V/J rows
326, 333, 369 and 748. No other original positive input is lost.

Every accepted IMGT/AHo map and chain matches the prior Immunum result. Martin
and Kabat recover the nine light-chain span failures described above, leaving
the unsupported terminal residue outside the numbered domain; all other accepted
maps remain identical. The 90 modest truncations retain parent consistency.
There are no unexplained numbering changes in this panel.

| Scheme | Exact maps / original valid fixture arrays | Identical labels / compared accepted residues |
| --- | ---: | ---: |
| IMGT | 1,454/1,483 | 169,566/169,645 |
| Martin | 1,456/1,484 | 169,509/169,645 |
| Kabat | 1,456/1,483 | 169,561/169,645 |
| AHo | 1,473/1,484 | 169,616/169,645 |

The exact-map denominators retain rejected inputs as nonmatches and keep the
original malformed-array exclusions. Lower totals than raw Immunum reflect
scope rejection, not changed accepted labels. These remain fixture-agreement
measures rather than independent accuracy estimates.

Every accepted input has a qualifying V match. Five lack the required J evidence:
numbering row 720 and COVID rows 4244, 4617, 10960 and 13044. They return
`j_match=None` with a diagnostic, preserving their numbering. No ancestral or
species-specific accuracy conclusion follows from these similarities.

Additional checks reconstruct the displayed query from the original input plus
imputed residues for all 1,569 accepted numbering/truncation rows in each scheme.
All 12,552 before/after renderings preserve sequence and the requested 80-column
width; repeating imputation preserves residue positions and amino acids.
Repository tests cover tied-reference disagreement, unavailable endpoints,
unknown symbols, partial references, mixed CDR definitions, multiple domains,
long-insertion conversion limits, numbered insertion ordering, immutable Python
results, and reference-dependent display order.

### Runtime and reuse

Measurements use Linux x86-64, an AMD Ryzen 9 9950X3D, Rust 1.96.0, CPython
3.13.13 and locked release builds. Each process was pinned to the same allowed
CPU. Imports and warm-up are excluded; timings are medians of seven serial
batches after five warm-up calls. No competing qualification process was running.
The parent is sequence-alignment PR #25 at `1ace91e`; feature timings use
`5b798a3`. These are workload-specific measurements, not cross-engine speed claims.

For the first 32 numbering-fixture inputs, align each sequence globally against
itself with input index 30 changed to Y (F when already Y). Each batch repeats
32 alignments 100 times; formatting repeats the first result 1,000 times at
width 80, plain color, with rulers.

| Operation | Parent | Antibody branch |
| --- | ---: | ---: |
| 32 pairwise sequence alignments | 1.818 ms | 1.807 ms |
| One pairwise alignment display | 3.08 μs | 7.51 μs |
| One numbered antibody display | — | 16.57 μs |

Calculation time is unchanged within measurement variability. The shared
annotated renderer adds about 4.4 μs to this pairwise display, reflecting its
per-cell coordinate labels and general row layout.

Numbering the same 32 inputs with eager V/J matching takes 295.58 ms per batch
(9.24 ms/input). An isolated build with only the private germline-matching call
omitted takes 4.46 ms (0.139 ms/input), retaining the recognition and conversion
work. This isolated measurement preceded the public [matching opt-out](#optional-germline-matching).
Matching accounts for about 98.5% of this workload's elapsed time.

For broader chain coverage, take the first ten H, K and L numbering fixtures
according to the prior Immunum results. Each batch repeats ten numbering calls
20 times, searching all reference species:

| Chain | Median per input |
| --- | ---: |
| H | 9.10 ms |
| K | 2.19 ms |
| L | 1.08 ms |

Bundled profiles and encoded references are reused across calls. Aligning ten
already-numbered heavy chains takes 69.65 μs; imputing the first stored result
takes 8.51 μs, each measured over 1,000 calls per batch. Neither operation repeats
V/J alignment. The full four-scheme quality run used concurrent processes and
therefore supplies no additional comparative timing result.

### Packaging and final checks

The final `a9f1337` source and parent were built on the same machine with
`maturin build --release --features python --locked` and
`cargo build --release --locked`, without an additional strip step. Sizes are
bytes; gzip measurements compress the executable alone at level 9 with zero mtime,
not a release archive containing extra files.

| Artifact | Parent | Antibody branch | Increase |
| --- | ---: | ---: | ---: |
| CPython 3.13 wheel | 9,606,369 | 9,893,436 | 287,067 (2.99%) |
| Python extension, unpacked | 34,904,672 | 36,107,808 | 1,203,136 (3.45%) |
| CLI executable | 51,288,144 | 52,561,568 | 1,273,424 (2.48%) |
| CLI executable, gzip | 12,333,380 | 12,616,583 | 283,203 (2.30%) |

The three-species reference subset at this measurement was 255,885 bytes. The lockfile adds Immunum 1.3.1,
strum 0.27.2 and strum_macros 0.27.2 without upgrading existing dependencies.
Immunum's default features are disabled; its CLI and Python layers are unused.

The installed wheel numbers antibodies, finds V/J matches, and constructs an
antibody alignment offline. It includes IMGT attribution; the source distribution
includes both attribution and the reference FASTA. Both exclude `docs/`; the
Rust package file list also retains the reference assets and excludes docs.
The source distribution does not bundle benchmark fixtures.
CLI release archives include the same attribution alongside the executable.

Validation at `a9f1337` passed 199 Rust library tests, 3 binary tests, 17 CLI tests,
8 doctests, 21 Python tests, Python type checking, and the repository's required
format/lint checks. The [user guide](https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md)
describes the supported scope and remaining limitations.

## Rat and rabbit reference expansion

The same pinned IMGT release now supplies rat (268 V / 13 J) and rabbit
(123 V / 20 J) references, including strain metadata and H/K/L for both species.
The bundle totals 1,603 V / 91 J records. All added J records contain the
expected W/F-G anchor: heavy-chain coverage ends at IMGT 128, light-chain
coverage at 127. Missing reference positions remain unavailable to imputation.
Synthetic V/junction/J cases test matching, species restrictions and terminal
imputation for both species, every chain class and all four numbering schemes.
These checks exercise reference integration; they do not establish biological
numbering or germline-assignment accuracy.

On 11 September 2026, the same locked release Python extension compared
`species=["human", "mouse", "alpaca"]` with the expanded default search.
Use the first ten accepted H, K and L numbering fixtures in input order,
classified by the initial IMGT adapter results. Each of seven batches numbers
those ten sequences ten times per setting, alternating setting order between
batches after five warm-up calls per setting. The process was pinned to CPU 0
on the runtime host described above, with no concurrent test or benchmark load.
Times include numbering, eager matching and Python result construction, excluding
imports, warm-up and rendering; report the median per input.

| Chain | Three species | Five species | Increase |
| --- | ---: | ---: | ---: |
| H | 9.138 ms | 10.863 ms | 18.9% |
| K | 2.198 ms | 4.156 ms | 89.1% |
| L | 1.088 ms | 1.433 ms | 31.8% |

Numbering and V reference sets are unchanged in these 30 cases; six J reference
sets change. Expanding the search can add ties or change a best match, so callers
needing a fixed search scope should specify species. An unrestricted search
pays for the additional candidates; the measured increase is approximately
1.73 ms for H, 1.96 ms for K and 0.35 ms for L.

The raw reference file grows from 255,885 to 349,651 bytes (+93,766).
Compressing each exact bundled file with Python `gzip.compress`, level 9 and
`mtime=0`, gives 35,752 and 51,355 bytes (+15,603). These are data-file sizes,
not measured wheel or executable growth. Source records and their ordering are
reproduced by the [preparation script](https://github.com/y1zhou/arpeggia/blob/master/data/germlines/prepare.py);
the [attribution file](https://github.com/y1zhou/arpeggia/blob/master/data/germlines/README.md)
records the expanded subset checksum.

## Llama reference supplement

On 15 September 2026, the same 30-input, seven-batch protocol and CPU 0 affinity
compared explicit `species=["human", "mouse", "alpaca", "rat", "rabbit"]` with
the six-species default on one locked release build. The supplement adds six
heavy-chain V and five J references; light-chain candidates are unchanged.

| Chain | Five species | Six species | Change |
| --- | ---: | ---: | ---: |
| H | 10.833 ms | 10.902 ms | +0.64% |
| K | 4.151 ms | 4.150 ms | −0.03% |
| L | 1.433 ms | 1.432 ms | −0.07% |

Numbering and V reference sets are unchanged in all 30 cases; one J reference
set changes. Synthetic llama V/junction/J inputs pass matching and terminal
imputation checks in all four schemes, without establishing biological accuracy.
The [supplement](https://github.com/y1zhou/arpeggia/blob/master/data/germlines/README.md#llama-protein-display-supplement)
adds 1,327 raw bytes. Concatenated reference data grow from 51,355 to 51,551 bytes
under the same gzip settings above; these are not package-size measurements.

## Terminal framework regression follow-up

A terminal sweep on 11 September 2026 used the first ten accepted H, K and L
numbering fixtures whose original IMGT maps start at 1 and end at 127/128.
Remove 1, 5, 9, 20, 22, 23, 24, 25 or 26 N-terminal residues and three C-terminal
residues, then compare each surviving residue's label and region against its
full parent in all four schemes: 1,080 calls.

Checking conversion windows alone left six inconsistent outputs from two lambda
parents (numbering rows 791 and 845, cut by 22 residues) in IMGT, Martin and
Kabat. The raw profile alignment had already shifted surviving framework
residues into CDR1 after losing coverage of Cys23. Requiring the profile span
to include IMGT 23–118 rejects these fragments and also protects every count-based
core conversion window in the pinned backend. All 520 accepted results preserve
parent numbering and regions; the remaining 560 produce explicit coverage errors.
Each scheme accepts 130 and rejects 140. This is a conservative supported-input
boundary, not proof that every rejected fragment is intrinsically unnumberable.
The saved accepted IMGT maps from the original qualification already span this
interval; the new sweep exercises a boundary absent from that panel.

## Optional germline matching

On 12 September 2026, a locked release Python extension compared the default
five-species search with `match_germlines=False` on the same 30 inputs and host
as the reference-expansion benchmark. Use the same seven batches, ten repeats
of each chain's ten inputs, five warm-up calls per setting and CPU 0 affinity.
Alternate enabled/disabled order between batches; no test or build ran alongside
the measurements. Timings include numbering and Python result construction,
excluding imports, warm-up and display.

| Chain | Matching enabled | Numbering only | Speedup |
| --- | ---: | ---: | ---: |
| H | 10.848 ms | 0.135 ms | 80.2× |
| K | 4.149 ms | 0.122 ms | 34.1× |
| L | 1.429 ms | 0.122 ms | 11.7× |

All 30 inputs retain identical residue labels, regions, input correspondence,
domain spans, chain assignments, confidence and numbering diagnostics. Only
search state and V/J results differ. These are warm, per-input measurements
for the sampled chains; they do not measure process startup or promise a fixed
speedup for other inputs. Matching remains enabled by default.

Numbering-only runs also repeated the pinned 26,365-input qualification panel
in all four schemes, comparing against the stored matching-enabled results.
Every scheme retains all 25,474 accepted position maps, chain assignments,
domain spans and numbering diagnostics, with the same 891 rejections. Search
state is false and both matches are absent throughout; matching-specific
missing-reference diagnostics are omitted as expected.

## Display validation

The 11 September 2026 display revision checked original-input rulers across gaps
and wraps, continuous V/J counts, blank imputed coordinates, reference-defined
CDR bands, yellow imputation highlights, tied names and Unicode width in Rust,
CLI and Python. It passed 202 library tests, three binary tests, 18 CLI tests,
eight doctests, 21 Python tests, type checking and pre-commit hooks.

Six representative AntPack full/truncated inputs were rendered before and after
imputation under all four schemes: 48 displays preserved supplied and imputed
sequence content, respected width, and retained idempotent imputation. This
revision did not repeat the full numbering benchmark.
