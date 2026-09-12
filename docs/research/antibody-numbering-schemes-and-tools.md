# Antibody numbering schemes and high-throughput tools: September 2026 assessment

**Research dates:** Initial survey 8 September 2026; Immunum and IMGT implementation follow-up 10 September 2026 (section 9). Other package release records retain the initial survey's cutoff.
**Scope:** Antibody variable-domain numbering, structural correspondence and insertion placement, performance evidence, and Python/Rust integration.
**Evidence:** Primary papers, official package records, source files, and maintainer issue discussions.
**Validation status:** Section 9 records Rust integration probes and reference-data
counts. The [three-engine comparison](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md)
covers 26,365 protein inputs, measuring fixture agreement and recognition
behavior rather than independent accuracy or comparative throughput.

## Findings

A numbering scheme defines positional equivalence; an engine assigns residues
to those positions. Martin and AHo address structural limitations in older
conventions, but no scheme or engine establishes universally correct structural
correspondence.[^martin][^aho][^evaluation2024]

Arpeggia selected Immunum's Rust core after the fixture comparison. AntPack,
ANARCII and RIOT provide complementary CPU, difficult-case and protein V/J
comparators. The evidence below does not establish a universal speed/accuracy
winner; reported performance belongs to the cited versions and workloads.

## 1. Schemes, engines, regions and correspondence

| Concept | What it specifies | Example |
|---|---|---|
| **Numbering scheme** | Standard positions, allowed insertions/deletions, and a notion of equivalent positions | IMGT, Kabat, Chothia, Martin, AHo |
| **Numbering engine** | How a sequence is mapped to a scheme | HMM alignment, consensus alignment, or a learned sequence-to-label model |
| **CDR definition** | Which positions are called CDR versus framework | A scheme-associated convention or an independently selected region definition |
| **Structural correspondence** | Which observed residues occupy equivalent three-dimensional roles | Curated or computed correspondence between crystallographic domains |

The original ANARCI paper distinguishes schemes and implements several of them. The newer structural evaluation documents that differences in scheme and CDR definitions can change which residues are included in a loop.[^anarci][^evaluation2024] AntPack documentation also explicitly supports numbering with one convention while using another for region assignment in some APIs.[^antpack-regions]

For a reproducible humanization or design pipeline, save **both** `numbering_scheme` and `cdr_definition`. Do not use an unlabeled field such as `cdr1` as though its boundaries were self-evident. Similarly, “Chothia,” “enhanced Chothia/Martin,” and “AbM CDR definition” should not be collapsed into one unnamed convention.

## 2. Numbering schemes

### 2.1 Kabat and Chothia are different starting points

Kabat's sequence-oriented convention and Chothia's structural interpretation do not always place gaps or delimit loops in the same way. The practical issue is not merely whether the displayed position numbers look familiar: a different insertion placement changes which residues are treated as corresponding across antibodies.[^martin][^evaluation2024]

Preserve a legacy dataset’s original scheme when reproducing its annotations or mutations; familiarity alone does not establish structural accuracy.

### 2.2 Martin / enhanced Chothia: explicit structural corrections

Abhinandan and Martin's 2008 paper proposed structurally motivated corrections to existing numbering, including insertion placement.[^martin]

Martin/enhanced Chothia is therefore worth supporting when the workflow depends on structural modeling or comparison to structural antibody literature. The original ANARCI scheme support includes enhanced Chothia, so adopting this scheme does not itself require adopting a newly released engine.[^anarci]

### 2.3 AHo: a structure-derived coordinate system

Honegger and Plückthun's AHo scheme was developed from structural alignment of immunoglobulin variable domains. Its gap placement was chosen to support spatial correspondence across domain classes rather than simply preserve a historical sequence numbering convention.[^aho]

This makes AHo a particularly relevant candidate for cross-antibody structural analysis. It does not imply that every unusual loop has a unique spatial alignment, or that every tool assigning AHo numbers infers that alignment correctly. Scheme design and engine correctness remain separate questions.

### 2.4 IMGT: a useful interoperability default

IMGT provides an extensively documented position system across immunoglobulin and related receptor variable domains. Its standard V-domain boundaries are FR1 1–26, CDR1 27–38, FR2 39–55, CDR2 56–65, FR3 66–104, CDR3 105–117, and FR4 118–128. Conserved landmarks include positions 23, 41, 104, and 118.[^imgt]

IMGT provides a documented interchange convention; structural and legacy workflows can still require alternate mappings.

Insertion positions must be interpreted with IMGT's ordering rules, especially around long CDR3s. **Do not assume that sorting position labels as ordinary strings reconstructs sequence order.** The safe implementation stores explicit input offsets and a scheme-aware ordering rather than deriving order from a displayed label.[^imgt]

### 2.5 What the 2024 reassessment adds

Zhu, Olson, and Magliery's 2024 study, *50 Years of Antibody Numbering Schemes*, compared statistical variation and structural behavior across common conventions. It found differences in CDR coverage and limitations in structural correspondence, including light-chain loop alignment issues. Its findings argue against assuming that any one familiar CDR convention perfectly captures all structurally or statistically important residues.[^evaluation2024]

## 3. Insertion-placement failure modes

There are at least three distinct failure modes:

**An engine error within the selected scheme.** A residue is assigned to an incorrect framework position, a domain boundary is misplaced, or an insertion is assigned to the wrong permitted location. Improving the engine can fix this without changing the scheme.

**A scheme limitation.** The scheme's prescribed correspondence is not the best structural match for the particular loop or chain class. A perfectly compliant engine can still produce a geometrically imperfect alignment.

**Genuine ambiguity.** Repeated residues, extreme loop lengths, truncation, or unusual structures permit multiple defensible correspondences. An engine may need to return a flag or alternative rather than one overconfident label.

The structural papers and newer engine evaluations document examples motivating these distinctions.[^martin][^aho][^evaluation2024][^anarcii-paper] The proposed distinction is useful because it determines the remedy: change the engine, change the convention, or preserve uncertainty.

A correctly identified conserved cysteine is a useful diagnostic, not proof that every intervening residue has been numbered correctly. Conversely, a deliberately mutated anchor should be flagged for review rather than automatically declared impossible. Sequence identity at an anchor and structural location of an anchor are different tests.

## 4. Engine comparison

| Tool | Method / integration | Best reason to evaluate it | Main reservation |
|---|---|---|---|
| **ANARCI** | Reference HMMs; Python-oriented ecosystem | Legacy reproducibility and comparison to established annotations | Performance and outputs depend on implementation, reference set, batching, and optional germline work |
| **AntPack** | Optimized native alignment routines with Python APIs | Strong CPU-throughput baseline for conventional antibody numbering | Hard-case accuracy is not uniform; current versus legacy licensing differs |
| **ANARCII** | Sequence-to-label language model; Python/PyTorch | Difficult sequences, unusual formats, and GPU batches | Model inference and deployment costs; not a universal fastest CPU tool |
| **Immunum** | Rust core, Python, Polars, JavaScript/WASM | Native Rust integration and vectorized data processing | Independent accuracy evidence and comparative benchmarks need strengthening |
| **RIOT** | Broader nucleotide/amino-acid immunoglobulin annotation | Combine numbering with germline and repertoire annotation | Different scope from pure numbering; restrictive use conditions must be checked |
| **AbNumber** | High-level Python wrapper with selectable backend support in inspected source | Convenient position objects, slicing, and region-aware workflows | Wrapper behavior and backend must be specified; not a new inference algorithm |

Sources for the engine descriptions are the original papers and inspected project documentation/source.[^anarci][^antpack-paper][^antpack-doc][^anarcii-repo][^immunum-source][^riot][^abnumber]

### 4.1 AntPack

The peer-reviewed 2024 Parkinson and Wang paper reports a benchmark of 3,492 antibody sequences from structural data, repeated five times on an Intel i7-13700K with SSD storage. AntPack completed the numbering in **under 0.5 seconds**, compared with **35–45 seconds for ANARCI** and **12–13 seconds for AbRSA** in that setup.[^antpack-paper]

This supports a substantial throughput improvement over that ANARCI configuration. It does not establish the same multiplier against every parallelized ANARCI installation, against optional V/J assignment turned off, or against later alternatives. Agreement with other engines is also not automatically a structure-based gold standard.

The official project describes optimized position-specific alignment routines and provides Python APIs for single-chain and paired-chain processing. TCR and antibody workflows should not be assumed to have identical runtime or support characteristics.[^antpack-doc]

Use AntPack as a CPU baseline, retaining failures and stratifying unusual loops, truncations, ambiguous residues and uncommon formats.

### 4.2 ANARCII

ANARCII's peer-reviewed paper was published on **21 May 2026**, with an August version-of-record date. It reports a sequence-to-label model and advantages on difficult sequence classes. In **28 VNAR structures**, AntPack misplaced Cys104 in **12**, while ANARCII placed that anchor correctly in all 28. This is a targeted structural result, not proof of perfect full-domain numbering.[^anarcii-paper]

The paper also reports approximately **90,000 sequences/minute** for its speed model and **70,000/minute** for its accuracy model on an **A100 GPU**. Its ANARCI comparator reached about **75,000/minute on 32 CPUs**, without V/J assignment. These differently provisioned measurements do not establish a universal wall-clock or cost winner.[^anarcii-paper]

Its difficult “no-truth” dataset deliberately enriches disagreements between ANARCI reference versions; rates from that set should not be extrapolated to ordinary complete therapeutic antibodies. Some CDR2/DE-region assignments remained ambiguous.[^anarcii-paper]

For deployment, the verified PyPI release is **2.0.8**, dated **30 June 2026**, with **Python ≥3.11** specified. The inspected repository license is **BSD 3-Clause**. The project provides a Python package and documented model-based numbering workflow.[^anarcii-release][^anarcii-license][^anarcii-repo]

ANARCII merits difficult-case and GPU-batch comparisons. Scores require calibration before probability interpretation; multidomain recovery requires explicit validation.

### 4.3 Immunum

Immunum implements a semi-global Needleman–Wunsch-style alignment against position-specific consensus scoring matrices. The core is Rust, with Python bindings, a Polars plugin, and JavaScript/WASM support. It recognizes antibody heavy, kappa, and lambda chains plus four TCR chain classes. The repository is MIT-licensed.[^immunum-source]

This is an attractive engineering combination for a service or data pipeline: native Rust processing can be used directly, while Python users can stay in a tabular workflow rather than repeatedly calling a Python function for each sequence. Those are integration advantages, not independently measured accuracy advantages.

The 10 September follow-up verified **1.3.1, released 2 September 2026**, on both PyPI and crates.io. Its published Rust artifact identifies commit `45bb70d34802cc592ebd86e685cc9f551885a2d6`. Chothia, Martin and AHo were added on 28 August; version 1.3.1 contains all five schemes, with alternate antibody schemes derived from internal IMGT numbering.[^immunum-release][^immunum-commit][^immunum-source]

Pin the deployed artifact and test each required scheme and chain class. Section 9 records integration boundaries that release availability alone does not resolve.

#### Comparative benchmark limits

The public evidence does not justify a blanket conclusion. Two maintainer issues are particularly relevant:

| Evidence | What it says | Implication |
|---|---|---|
| Issue #32, “Fix `antpack` parallelization benchmark” | Maintainers identify problematic scaling in the AntPack parallel comparison and propose fixes | The comparative speed chart is not a settled basis for declaring a universal winner |
| Issue #33, “Fix correctness benchmarks” | The current reference uses positions where ANARCI and AntPack agree; independent structural truth is requested | Agreement on a consensus subset is not independent structural accuracy |

Both issues were open in the retrieved records.[^immunum-speed-issue][^immunum-truth-issue]
Our fixture comparison supports Immunum as the first-choice Rust candidate,
subject to the qualification issues in section 9; it establishes neither
independent structural accuracy nor a comparative speed ranking.

### 4.4 RIOT: separate V/J alignment for protein numbering

RIOT aligns protein V and J segments separately, then transfers their germline-to-scheme mappings onto the input. This avoids treating the intervening CDR3 as one long insertion against a complete V–J profile. Its paper illustrates an ANARCI failure on 4ocr; that example does not establish an Immunum failure.[^riot]

The reported protein V/J assignment results, **96.94%/97.72%** on 1,274 therapeutic sequences, measure agreement with exhaustive Smith–Waterman/E-value assignments, not residue-numbering accuracy. Numbering comparisons found 23 IMGT conflicts and left some disagreements unresolved. Immunum was not a comparator; the paper supplies no dedicated alpaca/VHH numbering-accuracy result.[^riot]

The direct protein comparison is recorded in section 9. The inspected source
has IMGT/Martin/Chothia/Kabat but no AHo; its Rust component performs prefiltering
while the numbering pipeline remains Python. Native integration convenience
does not establish numbering quality.[^riot-source]

### 4.5 AbNumber: useful API, but disclose the backend

The inspected AbNumber source supports an optional ANARCII path, performs scheme conversion, and can convert to a legacy-compatible result. It falls back to ANARCI when germline assignment is requested. Thus “we used AbNumber” is not enough to identify the inference engine.[^abnumber]

The source also contains defensive handling of duplicate positions returned by ANARCI. That is a useful reminder to validate output invariants even when a mature wrapper is used. For large jobs, verify whether the calling pattern reuses an engine or repeatedly constructs one; convenient per-chain APIs can otherwise obscure startup costs.

## 5. Release and licensing audit at the cutoff

| Package / source | Verified state | Deployment consequence |
|---|---|---|
| **AntPack 0.4** | Non-yanked default release shown on PyPI; released 6 August 2025 | Versions 0.3.9 onward have academic/noncommercial licensing and license-key setup |
| **AntPack 0.5** | Released 14 April 2026 but **yanked**, reason “Bug fix” | Do not treat a newer-looking release number or README as a safe default |
| **AntPack 0.3.8.6.3** | GPL line; released 23 June 2026 | A lower version can be newer by calendar date; assess the exact artifact and its GPL obligations |
| **ANARCII 2.0.8** | Released 30 June 2026; Python ≥3.11; repository BSD 3-Clause | Pin model/package versions and account for the PyTorch environment |
| **Immunum 1.3.1** | PyPI and crates.io release 2 September 2026; Rust artifact matches `45bb70d…` | All five schemes are released; disable default features for native core use |

These records were checked against PyPI and repository sources, not inferred from version ordering alone.[^antpack-current][^antpack-yanked][^antpack-gpl][^anarcii-release][^anarcii-license][^immunum-release][^immunum-commit]

Licensing here is a report of the maintainers' published terms, not legal advice about a specific deployment. Review the actual package, weights, bundled references, and dependencies. In particular, **do not describe all AntPack versions as GPL or all current AntPack versions as commercially unrestricted**.

## 6. Python and Rust examples

These are interface examples grounded in the inspected documentation/source. The follow-up compiled Immunum's native core; the Python examples were not executed.

### 6.1 Immunum in Python

```python
from immunum import Annotator

annotator = Annotator(chains=["H", "K", "L"], scheme="imgt")
sequence = (
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGG"
    "VIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYC"
    "AREGTTGKPIGAFAHWGQGTLVTVSS"
)
result = annotator.number(sequence)
print(result.chain)
print(result.confidence)
print(result.numbering)
```

For tabular processing, the documented Polars plugin exposes `imp.number(...)` and `imp.segment(...)` as expressions.[^immunum-source]

### 6.2 Immunum in Rust

```rust
use immunum::{Annotator, Chain, Scheme};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let annotator = Annotator::new(
        &[Chain::IGH, Chain::IGK, Chain::IGL],
        Scheme::IMGT,
        None,
    )?;

    let sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
    let numbered = annotator.number(sequence)?;
    println!("{:?}", numbered.chain);
    Ok(())
}
```

The follow-up compiled this constructor/method pattern against `immunum = { version = "=1.3.1", default-features = false }`, keeping upstream CLI, Python, Polars and WASM integrations out of Arpeggia's dependency graph.[^immunum-manifest]

### 6.3 ANARCII in Python

```python
from anarcii import Anarcii

# Reuse one model instance for a batch rather than constructing one per chain.
model = Anarcii()
imgt_results = model.number([sequence])
martin_results = model.to_scheme("martin")
```

This numbering/conversion pattern is also used by the inspected AbNumber backend. Select and record model mode, device, batching, and version according to the deployed ANARCII documentation; the example intentionally does not imply a particular throughput.[^abnumber][^anarcii-repo]

## 7. Future structural workflows

Sequence-only Arpeggia behavior is specified in
[ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md).
A future structure workflow additionally needs a correspondence between input
offsets, domain identity, numbered positions and observed-structure residue IDs.
Chain IDs, author residue numbers and insertion codes remain distinct from
sequence offsets. Missing coordinates must remain explicit.

Whether to number a declared complete sequence or only observed residues needs
a separate decision: deleting unresolved residues can alter apparent loop
lengths, while declared sequence has no coordinates for fitting. Numbering-based
framework/CDR selections would permit fitting one region and evaluating another.

Independent-engine escalation should target detection failures, weak evidence,
unusual insertions or truncations, inconsistent labels/CDR boundaries, multiple
domains and structurally misplaced anchors. Engine agreement is evidence only
to the extent their errors are independent. Numbering, germline similarity,
humanness and developability require separate validation and runtime accounting.

Scheme conversion preserves a mapping; it does not realign the input against a
scheme-specific structural reference. This matters when comparing engines that
share an IMGT intermediate.[^immunum-source][^abnumber]

## 8. Independent validation criteria

A useful test suite should stratify heavy/kappa/lambda, conventional antibodies/VHH/unusual formats, species, domain completeness, CDR lengths, ambiguous residues, framework insertions, and multi-domain constructs. Include the actual kinds of de novo or heavily engineered sequences expected in this project rather than relying only on standard human therapeutic antibodies.

| Dimension | What to measure | Common misleading substitute |
|---|---|---|
| Domain detection | Correct number of domains and correct spans | At least one numbered region returned |
| Scheme assignment | Residue-level labels and insertion ordering | Percentage agreement only where old tools agree |
| CDR annotation | Both start and end boundaries under an explicit convention | A correct conserved cysteine alone |
| Structural consistency | Curated correspondence in well-resolved regions; ambiguity retained | Treating any predicted coordinate as ground truth |
| Failure handling | Rejections, warnings, malformed labels, silent misassignments | Reporting accuracy only on successful outputs |
| Throughput | Cold and warm wall time, batch size, cores, device, memory | Mixing single-core CPU, 32-core CPU, and A100 results |
| Reproducibility | Exact software/model/reference versions and optional tasks | Only the tool's name |

For reference labels, use structurally curated cases where a correspondence is defensible and explicitly mark unresolved regions. Keep all examples out of tuning and avoid filtering away cases merely because engines disagree. The maintainers' correctness issue in Immunum is a particularly clear example of why an agreement-defined reference can otherwise bias the comparison.[^immunum-truth-issue]

## 9. Arpeggia implementation findings

Accepted scope and API decisions from the completed design interview are recorded in
[ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md).
The selected backend is Immunum's Rust core with default features disabled;
Arpeggia supplies its own result objects, Python bindings and CLI rendering.

### 9.1 Native engine qualification

Immunum 1.3.1 is the smallest integration candidate found: its core compiled
without Python and already implements IMGT, Martin, AHo and Kabat. The published
crate is 91,762 bytes; generated consensus JSON totaled 554,835 bytes in the
probe build. Neither figure measures the eventual Arpeggia binary or wheel
increase.[^immunum-release][^immunum-manifest]

The annotator returns one best variable domain, inclusive zero-based input
bounds and positions in input order. Arpeggia uses half-open spans. Its
`confidence` is a normalized alignment heuristic, not a calibrated probability.
Profiles derive from human/mouse RepSeqIO 1.9 V/J sequences; there is no species
filter or germline assignment. Alpaca accuracy needs separate qualification.
Construction parses embedded profiles; an internal `RefCell` buffer prevents
sharing one annotator across parallel workers without synchronization. Reusing
instances per worker is a candidate for avoiding repeated construction.[^immunum-annotator][^immunum-profiles]

A standalone debug-build probe used the 122-residue IGH sequence in section 6:

| Probe | Observed result |
|---|---|
| Number the sequence under IMGT, Martin, AHo and Kabat | All four returned results; this checks integration, not positional accuracy |
| Number `sequence + GGGGSGGGGS + sequence` | All four returned only the second domain, inclusive span `132..=253`, without a multidomain warning |
| Apply the IMGT CDR3 conversion rule to length 65 | 52 insertions with valid single-letter labels |
| Apply that rule to length 66 | Invalid label `112[` |
| Apply that rule to length 413 | Integer-overflow panic |

The last three checks call `number_with_rules` directly with the
`Insertion::Symmetric { left: 111, right: 112 }` rule from `IMGT_RULES`.
These are synthetic conversion tests, not end-to-end numbering of biological
long-CDR antibodies. Single-character insertion labels and unchecked `u8`
arithmetic require a supported boundary or correction **before** conversion;
postprocessing malformed output would not prevent the panic.[^immunum-numbering]

The upstream AHo CDR table is explicitly marked unverified and its CDR3 start
disagrees with its cited alternatives. Martin uses AbM region boundaries.
Numbering support therefore does not settle the CDR coloring policy.[^immunum-aho][^immunum-numbering]

The [AntPack-fixture comparison](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md)
tested Immunum alongside RIOT and AntPack. It found strongest fixture-label
agreement in Immunum, but also a light-chain Martin/Kabat span defect and
non-antibody inputs passing the default confidence threshold. The Arpeggia adapter adds recognition, span and conversion guards; its
[qualification](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md#arpeggia-adapter-qualification)
distinguishes these corrections from fixture agreement. A bespoke engine would
add profile, detection and insertion-rule maintenance without accuracy evidence.
Constant-region numbering remains separate work: IMGT defines a distinct
C-domain system, and this engine models variable domains only.[^imgt-constant]

Source inspection identifies a public integration path: load the H/K/L
`ScoringMatrix` profiles, call `align()` with reusable `AlignBuffer` storage,
retain the winning raw `Alignment`, then call `apply_numbering()`. This exposes
domain evidence and source-position runs before conversion. `Annotator::number()`
hides that alignment and converts before checking confidence, so validating
only its final result cannot prevent the insertion failure. The lower-level
path reuses upstream scoring, alignment and scheme rules; it needs only the
chain-selection and validation orchestration.[^immunum-core-api]

The adapter's profile reuse, confidence/evidence checks, per-rule insertion
limits and multidomain rejection are recorded in
[ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md#inputs-and-result-objects).
These guards run before conversion to prevent upstream overflow or malformed labels.

Conversion boundaries need separate treatment from CDR membership. Martin/Kabat
can omit a terminal light-chain source position; the high-level AHo path can
append light-chain position 149 from the next input residue. Arpeggia preserves
the unnumbered input suffix and avoids duplicating an existing AHo 149.
Raw-alignment conversion also permits a separate CDR convention without
realignment.[^immunum-annotator][^immunum-numbering]

### 9.2 Germline data and interpretation

The IMGT terms retrieved on 10 September 2026 license data and metadata under
CC BY 4.0; tools retain separate terms. The attributed source snapshot avoids
reusing transformed datasets with older terms.[^imgt-terms]

Downloaded IMGT/GENE-DB release **202636-7** on 10 September 2026; the release and
amino-acid files reported last modification on 5 September. The gapped
`IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP` file was
3,332,989 bytes, SHA-256
`3cb6b0b8cb8940b3b2a9b105771a6a74aa67c06e3ca39eaea0d2030c90e7efd0`.[^imgt-download]

The bundled subset retains IGHV/IGKV/IGLV and IGHJ/IGKJ/IGLJ, functional
records including bracketed/parenthesized `F`, and species names with their
strain/subspecies suffixes. It excludes stop-containing sequences and retains
partial records, original IMGT gaps and duplicate source records:

| Species | V references | J references | Total |
|---|---:|---:|---:|
| Human | 511 | 33 | 544 |
| Mouse, including strains/subspecies | 628 | 19 | 647 |
| Alpaca | 73 | 6 | 79 |
| Rat | 268 | 13 | 281 |
| Rabbit | 123 | 20 | 143 |

The bundled five-species subset contains 1,603 V and 91 J records and occupies
349,651 bytes raw. Its release, checksum, attribution and exact transformation
are recorded with the [runtime data](https://github.com/y1zhou/arpeggia/blob/master/data/germlines/README.md).
Strain/subspecies suffixes are retained; exact species-name matching would
discard most mouse and rat references. The
[expansion benchmark](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md#rat-and-rabbit-reference-expansion)
compares three- and five-species search costs on the same build.

In the original human/mouse/alpaca subset of 1,212 V and 58 J records, all retained J
references lack partial-record flags and are 12–20 amino acids long; the alpaca
IGHJ5*01 sequence nevertheless ends before IMGT 128. FR4 coverage therefore
follows the observed W/F118 anchor and available residues, not a completeness
flag or assumed suffix length. V references span 48–105
ungapped residues; 45 have missing 5′ sequence represented by leading padding.
Four human V references contain six `X` characters. Preserve reference coverage
separately from internal alignment gaps, and exclude ambiguous residues from
known-residue evidence and imputation. Ordinary `SeqAlignment` identity retains
its existing literal-symbol semantics.

The snapshot contains alpaca heavy-chain V/J references but no alpaca K/L
references. Of the 73 retained alpaca V references, 66 are marked partial at the
3′ end around Cys104; discarding all partial entries would remove most coverage.
This supports an alpaca **VHH** germline scope, not exhaustive conventional
alpaca light-chain matching.[^imgt-alpaca]

Prefer separate V and J similarity results, retaining tied gene/allele names
and measured coverage. A stitched V+J display is not an inferred ancestral
antibody: junctional additions and D contributions are not recovered by that
operation. NCBI's IgBLAST reports D/J assignment only for nucleotide searches;
Arpeggia therefore labels its amino-acid J comparison as similarity, not a
unique gene call. Missing reference coverage must remain distinguishable from
a true deletion.[^igblast]

The accepted matching policy uses local BLOSUM62 with gap costs 10/0.5 in
separate V and J windows, with known-residue coverage gates. This reuses
Arpeggia's alignment contract without introducing uncalibrated E-values.
See [ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md)
for the windows, thresholds and missing-match behavior.

Protein sequences cannot distinguish some alleles: for example, retained human
IGHJ4*01/*02/*03 have identical amino-acid sequences. Reference identity and
provenance must survive sequence deduplication. Selecting a display representative
must not discard tied references used for imputation.

### 9.3 Position correspondence and display reuse

At `cbae1717f8bdb87fd504066dd6e3e3730de0dfe9`, antid delegates numbering, ordered
positions and V/J matching to AntPack. Its useful alignment model is the union
of numbered positions with one gapped sequence per row, not a new multiple
sequence alignment calculation. Preserve input order and validate compatible
schemes/chain classes. Its displayed germline merges V and J, favoring J where
both cover a position; that convenience should not conceal junction provenance.[^antid-numbering]

Rust can retain sequence bytes and typed position/input-offset records, deriving
region slices and tabular views when requested. Avoid duplicated per-residue
DataFrames. Position order must follow the scheme: for example, IMGT's inserted
positions on the 112 side run in reverse insertion order.[^imgt]

Arpeggia's `SeqAlignment` promises an optimal pairwise BLOSUM62 alignment.
Numbering-based correspondence need not be that optimum. Reuse the sequence
module's width handling, escaped names, rulers, foreground operation colors
and upstream `Style` backgrounds through a small internal renderer. Do not
construct a dummy `SeqAlignment` with misleading scores. Actual V/J pairwise
comparisons can still call `align_seqs()`. CDR backgrounds require named region
boundaries and handling of gap/blank cells. Region identification must remain
possible without color.[^arpeggia-alignment]

The implemented [display contract](https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md#display-and-antibody-alignments)
keeps source coverage and provenance separate from the stitched display row.

### 9.4 Comparing engines and CDR definitions

No inspected evidence establishes a RIOT-versus-Immunum accuracy ranking.
Immunum's maintainers identify agreement between ANARCI and AntPack as their
current benchmark reference and request independent structural labels.
Its 1.3.1 tests include an ultralong bovine CDR3 (4k3e H), but assert domain
coverage rather than every numbered position. ANARCI's limitations cannot be
assumed to apply unchanged to Immunum's position-dependent scoring.[^immunum-truth-issue][^immunum-alignment]

The [three-engine benchmark](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md)
separates returned mappings, residue-label agreement and domain coverage on
26,365 AntPack-fixture inputs. Its labels do not resolve independent accuracy;
that comparison did not test germline matching or species-specific structural
accuracy. Adapter and reference-integration checks are recorded separately.
Assess remaining disagreements against scheme rules and structural evidence.

A numbering scheme labels residues; a CDR definition assigns region membership.
For example, under fixed Chothia/Martin numbering, Kabat CDR-H1 is H31–H35 and
AbM/Martin CDR-H1 is H26–H35. H28 keeps its number but changes region. Mixed
definitions must be transferred through residue correspondence, not applied as
numeric cutoffs in another scheme.[^martin-chapter]

The accepted defaults are recorded in
[ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md).
The Martin group's 2024 study supports pairing AbM boundaries with Martin
numbering. The selected AHo structural-loop convention agrees with the core
boundaries illustrated by Honegger and Plückthun, but is not universal across
tools and differs from Immunum's unverified table. Public docstrings cite
these sources.[^martin-loops][^aho-loops][^aho-original]

Immunum's explicit Chothia CDR table follows the 2021 consensus: heavy
26–32 / 52–56 / 96–101 and light 26–32 / 50–52 / 91–96. These differ from its
Martin/AbM boundaries. Arpeggia exposes the distinct definition through explicit
`cdr_definition="chothia"`; the numbering argument's `chothia` alias still
resolves to Martin. This preserves intentional mixed conventions.[^immunum-chothia]

### 9.5 Explicit germline imputation

Antid's `imputed_seq` fills FR1/FR4 gaps from the first displayed germline and
requires the original variable-domain sequence to remain contiguous. It returns
a string, does not replace `X`, and does not fill CDRs.[^antid-imputation]

An unoccupied numbered position is not necessarily missing input. An eight-residue
IMGT CDR1 normally leaves positions 31–34 empty; filling them changes its loop
length. With raw amino-acid input alone, distinguishing an internal biological
deletion from omitted data requires information the string does not carry.[^imgt]

These distinctions motivate [terminal-only, provenance-preserving imputation](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md#terminal-imputation).

## Arpeggia integration

The implementation uses Immunum's raw alignment and guarded conversion APIs,
shared sequence rendering, and bundled IMGT references. Accepted behavior is in
[ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md);
arguments and display controls are in the
[user guide](https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md).
The [qualification report](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md#arpeggia-adapter-qualification)
retains adapter checks, display validation, runtime and package measurements.

## References and implementation records

[^martin]: Abhinandan KR, Martin ACR. **Analysis and improvements to Kabat and structurally correct numbering of antibody variable domains.** *Molecular Immunology* 45, 3832–3839 (2008). DOI: <https://doi.org/10.1016/j.molimm.2008.05.022>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/18614234/>.
[^aho]: Honegger A, Plückthun A. **Yet another numbering scheme for immunoglobulin variable domains: an automatic modeling and analysis tool.** *Journal of Molecular Biology* 309, 657–670 (2001). DOI: <https://doi.org/10.1006/jmbi.2001.4662>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/11397087/>.
[^evaluation2024]: Zhu Z, Olson KS, Magliery TJ. **50 Years of Antibody Numbering Schemes: A Statistical and Structural Evaluation Reveals Key Differences and Limitations.** *Antibodies* 13, 99 (2024). DOI: <https://doi.org/10.3390/antib13040099>. Publisher: <https://www.mdpi.com/2073-4468/13/4/99>; PubMed: <https://pubmed.ncbi.nlm.nih.gov/39727482/>.
[^anarci]: Dunbar J, Deane CM. **ANARCI: antigen receptor numbering and receptor classification.** *Bioinformatics* 32, 298–300 (2016; online 2015). DOI: <https://doi.org/10.1093/bioinformatics/btv552>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/26424857/>.
[^antpack-regions]: AntPack, official **clustering / region assignment** documentation, including separate numbering and CDR conventions: <https://antpackdocumentationlatest.pages.dev/clustering_overview>.
[^imgt]: IMGT, **IMGT unique numbering for V domains**, official scientific chart: <https://imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html>.
[^anarcii-paper]: Greenshields-Watson A, Agarwal P, Robinson SA et al. **ANARCII enables alignment-free antigen receptor numbering using a generalised language model.** *Communications Biology* 9, 1085 (2026). Published 21 May 2026; version-of-record date 12 August 2026. DOI: <https://doi.org/10.1038/s42003-026-10186-z>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/42162238/>. Numerical comparisons in this report are author-reported, not reproduced here.
[^antpack-paper]: Parkinson J, Wang W. **For antibody sequence generative modeling, mixture models may be all you need.** *Bioinformatics* 40, btae278 (2024). DOI: <https://doi.org/10.1093/bioinformatics/btae278>. This paper includes the AntPack numbering benchmark: <https://academic.oup.com/bioinformatics/article/40/5/btae278/7656770>.
[^antpack-doc]: AntPack, official **numbering background** documentation: <https://antpackdocumentationlatest.pages.dev/numbering_background>. The documentation site's displayed version can lag package releases; verify the deployed API separately.
[^anarcii-repo]: ANARCII, inspected release README and official user guide: <https://github.com/oxpig/ANARCII/blob/e0d8f192f5a861e03a50918f114d0f5735e42333/README.md>; <https://github.com/oxpig/ANARCII/wiki>.
[^immunum-source]: ENPICOM, **Immunum 1.3.1 README**, including algorithm, interfaces, chain/scheme support, and MIT license declaration: <https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/README.md>.
[^riot]: Dudzic et al. **RIOT—Rapid Immunoglobulin Overview Tool—annotation of nucleotide and amino acid immunoglobulin sequences using an open germline database.** *Briefings in Bioinformatics* 26, bbae632 (2025 issue; online 2024). DOI: <https://doi.org/10.1093/bib/bbae632>. Publisher: <https://academic.oup.com/bib/article/26/1/bbae632/7914577>; official project: <https://github.com/NaturalAntibody/riot_na>.
[^abnumber]: AbNumber, inspected **`abnumber/common.py`**, `_anarci_align` backend selection, scheme conversion, germline fallback, and duplicate-position handling: <https://github.com/prihoda/AbNumber/blob/master/abnumber/common.py>. File blob hash at retrieval: `d494cb3593745ee824727537176da5eec67e1aaf`. This is an inspected source record, not a guarantee about every packaged version.
[^anarcii-release]: ANARCII **2.0.8** package metadata and release history: <https://pypi.org/project/anarcii/>. Release date 30 June 2026; Python ≥3.11 in the retrieved package metadata.
[^anarcii-license]: ANARCII, pinned **BSD 3-Clause** license: <https://github.com/oxpig/ANARCII/blob/e0d8f192f5a861e03a50918f114d0f5735e42333/LICENCE>.
[^immunum-release]: Immunum **1.3.1**, released 2 September 2026: <https://pypi.org/project/immunum/1.3.1/>; Rust release metadata: <https://crates.io/api/v1/crates/immunum>. The inspected Rust artifact's VCS record identifies `45bb70d34802cc592ebd86e685cc9f551885a2d6`.
[^immunum-commit]: Immunum commit **`e027d5fe2405300508eee7ce78582a2fe4990f20`**, 28 August 2026, adding Chothia, Martin, and AHo and changing source version to 1.3.0: <https://github.com/ENPICOM/immunum/commit/e027d5fe2405300508eee7ce78582a2fe4990f20>.
[^immunum-speed-issue]: Immunum maintainer issue **#32**, *Fix `antpack` parallelization benchmark*: <https://github.com/ENPICOM/immunum/issues/32>. Opened 23 March 2026, updated 8 April 2026; open in the retrieved record.
[^immunum-truth-issue]: Immunum maintainer issue **#33**, *Fix correctness benchmarks*: <https://github.com/ENPICOM/immunum/issues/33>. Opened 23 March 2026; open in the retrieved record. The issue explicitly requests an independent, for example structure-based, gold standard.
[^riot-source]: RIOT [scheme and species enums](https://github.com/NaturalAntibody/riot_na/blob/2ee4dc3dcfa440cf04356d2c89dbf4917194d8dc/riot_na/data/model.py), [numbering pipeline](https://github.com/NaturalAntibody/riot_na/blob/2ee4dc3dcfa440cf04356d2c89dbf4917194d8dc/riot_na/api/riot_numbering.py), and [Rust manifest](https://github.com/NaturalAntibody/riot_na/blob/2ee4dc3dcfa440cf04356d2c89dbf4917194d8dc/Cargo.toml).
[^antpack-current]: AntPack, current PyPI package metadata and licensing: <https://pypi.org/project/antpack/>. At retrieval, the default non-yanked release was 0.4; versions 0.3.9 onward use academic/noncommercial terms and key setup.
[^antpack-yanked]: AntPack **0.5** release record and release history: <https://pypi.org/project/antpack/0.5/>. Released 14 April 2026; yanked for “Bug fix.”
[^antpack-gpl]: AntPack **0.3.8.6.3**, release record and GPL metadata: <https://pypi.org/project/antpack/0.3.8.6.3/>. Released 23 June 2026; Python ≥3.8 stated in this package record.
[^immunum-manifest]: Immunum 1.3.1 **Cargo manifest**: <https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/Cargo.toml>.
[^immunum-annotator]: Immunum **annotator**, including bounds, thresholds and scratch storage: <https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/annotator.rs>.
[^immunum-profiles]: Immunum **consensus profile provenance**: <https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/resources/consensus/README.md>.
[^immunum-numbering]: Immunum **numbering rules and insertion generation**: <https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/numbering.rs>.
[^immunum-aho]: Immunum **AHo numbering and unverified region table**: <https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/numbering/aho.rs>.
[^imgt-constant]: IMGT **unique numbering for C domains**: <https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVCsuperfamily.html>.
[^immunum-core-api]: Immunum 1.3.1 public [alignment state and API](https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/alignment.rs#L35) and [scheme conversion](https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/numbering.rs#L79).
[^imgt-terms]: IMGT **terms of use**, retrieved 10 September 2026: <https://www.imgt.org/about/termsofuse.php>.
[^imgt-download]: IMGT **GENE-DB downloads**: <https://www.imgt.org/download/GENE-DB/>; [release](https://www.imgt.org/download/GENE-DB/RELEASE); [gapped amino-acid references](https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F%2BORF%2BinframeP).
[^imgt-alpaca]: IMGT direct alpaca exports: [IGHV](https://www.imgt.org/genedb/GENElect?query=7.3+IGHV&species=Vicugna+pacos), [IGHJ](https://www.imgt.org/genedb/GENElect?query=7.6+IGHJ&species=Vicugna+pacos). Unfiltered exports contain 84 V and 7 J records; section 9 counts use the stated functional filter.
[^igblast]: NCBI **IgBLAST introduction**, including separate amino-acid and nucleotide capabilities: <https://www.ncbi.nlm.nih.gov/igblast/intro.html>.
[^antid-numbering]: antid **numbering objects, germline display and alignment**, inspected at `cbae1717f8bdb87fd504066dd6e3e3730de0dfe9`: <https://github.com/y1zhou/antid/blob/cbae1717f8bdb87fd504066dd6e3e3730de0dfe9/src/antid/numbering/antibody.py>.
[^arpeggia-alignment]: Arpeggia [sequence-alignment contract](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0009-separate-sequence-correspondence-from-rmsd-evaluation.md) and [renderer](https://github.com/y1zhou/arpeggia/blob/master/src/seq_alignment/display.rs).
[^immunum-alignment]: Immunum 1.3.1 [alignment implementation and 4k3e H coverage regression](https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/alignment.rs#L486).
[^martin-chapter]: Martin ACR, **Protein sequence and structure analysis of antibody variable domains**, Table 3.4: [author's chapter](https://citeseerx.ist.psu.edu/document?doi=ad292f8ef5c09a4ebb540a69e2c1e91366b1a3f2&repid=rep1&type=pdf).
[^martin-loops]: Martin group, **Do antibody CDR loops change conformation upon binding?** (2024), Martin numbering with AbM CDR boundaries: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10939163/>.
[^aho-loops]: Xu et al., **Functional clustering of B cell receptors using sequence and structural features** (2019), Methods: BCR notation: <https://pubs.rsc.org/en/content/articlehtml/2019/me/c9me00021f>.
[^aho-original]: Honegger and Plückthun (2001), original AHo paper, Figure 3: [author-hosted PDF](https://plueckthun.bioc.uzh.ch/wp-content/uploads/Publications/APpub0204.pdf).
[^immunum-chothia]: Immunum 1.3.1 [Chothia CDR definitions and cited 2021 consensus](https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/numbering/chothia.rs#L16).
[^antid-imputation]: antid [terminal imputation](https://github.com/y1zhou/antid/blob/cbae1717f8bdb87fd504066dd6e3e3730de0dfe9/src/antid/numbering/antibody.py#L592-L642).
