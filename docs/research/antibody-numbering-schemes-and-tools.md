# Antibody numbering schemes and high-throughput tools: September 2026 assessment

**Research cutoff:** 8 September 2026; this report does not assume developments later in September.
**Scope:** Antibody variable-domain numbering, structural correspondence and insertion placement, performance evidence, and Python/Rust integration.
**Evidence:** Primary papers, official package records, source files, and maintainer issue discussions.
**Validation status:** No numbering-engine accuracy or throughput benchmark was executed in this session. Published and maintainer-reported measurements are identified as such.

## Executive assessment

There have been meaningful improvements, but two different questions need separate answers. **A numbering scheme defines positional equivalence; a numbering engine assigns a sequence to those positions.** Replacing ANARCI with a faster engine does not automatically improve the underlying scheme, and converting an incorrect initial mapping into a different scheme does not necessarily correct it.[^anarci][^immunum-source]

For structural correspondence, **Martin/enhanced Chothia and AHo are important improvements over simpler historical conventions**, although they are not new in 2026. A 2024 structural/statistical study still found limitations among commonly used schemes rather than establishing one universally correct successor.[^martin][^aho][^evaluation2024]

For engines, my current shortlist is **AntPack for fast CPU numbering of conventional antibody sequences**, **ANARCII for difficult or unusual sequences and GPU batch processing**, and **Immunum for a Rust-native/Python/Polars deployment**. RIOT merits attention when numbering is part of a broader nucleotide/amino-acid germline-annotation workflow.[^antpack-paper][^anarcii-paper][^immunum-source][^riot]

**I did not find adequate evidence to declare that a newer tool universally surpasses AntPack in both speed and accuracy.** ANARCII provides published evidence of advantages on challenging sequence classes. Immunum is a serious engineering alternative, but its maintainers explicitly identify weaknesses in their AntPack parallel benchmark and in the independence of their correctness benchmark.[^anarcii-paper][^immunum-speed-issue][^immunum-truth-issue]

A deployment decision must also account for licensing and version differences: current AntPack licensing is not the same as its older GPL line, and Immunum's five-scheme source support is newer than the verified PyPI package.[^antpack-current][^antpack-gpl][^immunum-release][^immunum-commit]

## 1. Four concepts that should not be conflated

| Concept | What it specifies | Example |
|---|---|---|
| **Numbering scheme** | Standard positions, allowed insertions/deletions, and a notion of equivalent positions | IMGT, Kabat, Chothia, Martin, AHo |
| **Numbering engine** | How a sequence is mapped to a scheme | HMM alignment, consensus alignment, or a learned sequence-to-label model |
| **CDR definition** | Which positions are called CDR versus framework | A scheme-associated convention or an independently selected region definition |
| **Structural correspondence** | Which observed residues occupy equivalent three-dimensional roles | Curated or computed correspondence between crystallographic domains |

The original ANARCI paper distinguishes schemes and implements several of them. The newer structural evaluation documents that differences in scheme and CDR definitions can change which residues are included in a loop.[^anarci][^evaluation2024] AntPack documentation also explicitly supports numbering with one convention while using another for region assignment in some APIs.[^antpack-regions]

For a reproducible humanization or design pipeline, save **both** `numbering_scheme` and `cdr_definition`. Do not use an unlabeled field such as `cdr1` as though its boundaries were self-evident. Similarly, “Chothia,” “enhanced Chothia/Martin,” and “AbM CDR definition” should not be collapsed into one unnamed convention.

## 2. Have the established schemes improved structurally?

### 2.1 Kabat and Chothia are different starting points

Kabat's sequence-oriented convention and Chothia's structural interpretation do not always place gaps or delimit loops in the same way. The practical issue is not merely whether the displayed position numbers look familiar: a different insertion placement changes which residues are treated as corresponding across antibodies.[^martin][^evaluation2024]

For legacy datasets, maintaining the originally used scheme can be essential for reproducing annotations or mutations. For a new computational system, legacy popularity alone is not a sufficient reason to treat that scheme as a structural gold standard. That is a workflow recommendation, not a claim that one historical system has become invalid.

### 2.2 Martin / enhanced Chothia: explicit structural corrections

Abhinandan and Martin's 2008 paper directly addressed structural problems in existing numbering, including insertion placement, and proposed improvements. This is a substantive answer to the question about schemes that better fit antibody crystal structures: **structurally motivated correction of numbering was already an explicit research objective, not an innovation first introduced by the newest high-throughput tools.**[^martin]

Martin/enhanced Chothia is therefore worth supporting when the workflow depends on structural modeling or comparison to structural antibody literature. The original ANARCI scheme support includes enhanced Chothia, so adopting this scheme does not itself require adopting a newly released engine.[^anarci]

### 2.3 AHo: a structure-derived coordinate system

Honegger and Plückthun's AHo scheme was developed from structural alignment of immunoglobulin variable domains. Its gap placement was chosen to support spatial correspondence across domain classes rather than simply preserve a historical sequence numbering convention.[^aho]

This makes AHo a particularly relevant candidate for cross-antibody structural analysis. It does not imply that every unusual loop has a unique spatial alignment, or that every tool assigning AHo numbers infers that alignment correctly. Scheme design and engine correctness remain separate questions.

### 2.4 IMGT: a useful interoperability default

IMGT provides an extensively documented position system across immunoglobulin and related receptor variable domains. Its standard V-domain boundaries are FR1 1–26, CDR1 27–38, FR2 39–55, CDR2 56–65, FR3 66–104, CDR3 105–117, and FR4 118–128. Conserved landmarks include positions 23, 41, 104, and 118.[^imgt]

I would use IMGT as the canonical sequence-level interchange scheme in a new multi-tool pipeline, while retaining alternate scheme mappings where they serve a specific structural or legacy purpose. That recommendation prioritizes documented interoperability, not a claim that IMGT is the best geometric correspondence for every loop.

Insertion positions must be interpreted with IMGT's ordering rules, especially around long CDR3s. **Do not assume that sorting position labels as ordinary strings reconstructs sequence order.** The safe implementation stores explicit input offsets and a scheme-aware ordering rather than deriving order from a displayed label.[^imgt]

### 2.5 What the 2024 reassessment adds

Zhu, Olson, and Magliery's 2024 study, *50 Years of Antibody Numbering Schemes*, compared statistical variation and structural behavior across common conventions. It found differences in CDR coverage and limitations in structural correspondence, including light-chain loop alignment issues. Its findings argue against assuming that any one familiar CDR convention perfectly captures all structurally or statistically important residues.[^evaluation2024]

I did not identify a broadly adopted 2025–2026 scheme that resolves all these issues and replaces IMGT, AHo, and Martin. The defensible conclusion is **continued need for explicit conventions and structural validation**, not that recent tooling has made the scheme problem disappear.

## 3. What “correct insertion placement” actually requires

There are at least three distinct failure modes:

**An engine error within the selected scheme.** A residue is assigned to an incorrect framework position, a domain boundary is misplaced, or an insertion is assigned to the wrong permitted location. Improving the engine can fix this without changing the scheme.

**A scheme limitation.** The scheme's prescribed correspondence is not the best structural match for the particular loop or chain class. A perfectly compliant engine can still produce a geometrically imperfect alignment.

**Genuine ambiguity.** Repeated residues, extreme loop lengths, truncation, or unusual structures permit multiple defensible correspondences. An engine may need to return a flag or alternative rather than one overconfident label.

The structural papers and newer engine evaluations document examples motivating these distinctions.[^martin][^aho][^evaluation2024][^anarcii-paper] The proposed distinction is useful because it determines the remedy: change the engine, change the convention, or preserve uncertainty.

A correctly identified conserved cysteine is a useful diagnostic, not proof that every intervening residue has been numbered correctly. Conversely, a deliberately mutated anchor should be flagged for review rather than automatically declared impossible. Sequence identity at an anchor and structural location of an anchor are different tests.

## 4. Engine shortlist and what each option is actually good for

| Tool | Method / integration | Best reason to evaluate it | Main reservation |
|---|---|---|---|
| **ANARCI** | Reference HMMs; Python-oriented ecosystem | Legacy reproducibility and comparison to established annotations | Performance and outputs depend on implementation, reference set, batching, and optional germline work |
| **AntPack** | Optimized native alignment routines with Python APIs | Strong CPU-throughput baseline for conventional antibody numbering | Hard-case accuracy is not uniform; current versus legacy licensing differs |
| **ANARCII** | Sequence-to-label language model; Python/PyTorch | Difficult sequences, unusual formats, and GPU batches | Model inference and deployment costs; not a universal fastest CPU tool |
| **Immunum** | Rust core, Python, Polars, JavaScript/WASM | Native Rust integration and vectorized data processing | Independent accuracy evidence and comparative benchmarks need strengthening |
| **RIOT** | Broader nucleotide/amino-acid immunoglobulin annotation | Combine numbering with germline and repertoire annotation | Different scope from pure numbering; restrictive use conditions must be checked |
| **AbNumber** | High-level Python wrapper with selectable backend support in inspected source | Convenient position objects, slicing, and region-aware workflows | Wrapper behavior and backend must be specified; not a new inference algorithm |

Sources for the engine descriptions are the original papers and inspected project documentation/source.[^anarci][^antpack-paper][^antpack-doc][^anarcii-repo][^immunum-source][^riot][^abnumber]

### 4.1 AntPack: still an important performance baseline

The peer-reviewed 2024 Parkinson and Wang paper reports a benchmark of 3,492 antibody sequences from structural data, repeated five times on an Intel i7-13700K with SSD storage. AntPack completed the numbering in **under 0.5 seconds**, compared with **35–45 seconds for ANARCI** and **12–13 seconds for AbRSA** in that setup.[^antpack-paper]

This supports a substantial throughput improvement over that ANARCI configuration. It does not establish the same multiplier against every parallelized ANARCI installation, against optional V/J assignment turned off, or against later alternatives. Agreement with other engines is also not automatically a structure-based gold standard.

The official project describes optimized position-specific alignment routines and provides Python APIs for single-chain and paired-chain processing. TCR and antibody workflows should not be assumed to have identical runtime or support characteristics.[^antpack-doc]

My recommendation is to keep AntPack as a strong CPU baseline on representative conventional antibody data. Retain failure metadata and evaluate unusually long loops, truncations, ambiguous residues, and uncommon receptor formats separately. A fast successful return is not the same as a fully correct mapping.

### 4.2 ANARCII: the clearest recent published accuracy advance

ANARCII's peer-reviewed paper was published on **21 May 2026**, with an August version-of-record date. It reports a sequence-to-label model and advantages on difficult sequence classes. In **28 VNAR structures**, AntPack misplaced Cys104 in **12**, while ANARCII placed that anchor correctly in all 28. This is a targeted structural result, not proof of perfect full-domain numbering.[^anarcii-paper]

The paper also reports approximately **90,000 sequences/minute** for its speed model and **70,000/minute** for its accuracy model on an **A100 GPU**. Its ANARCI comparator reached about **75,000/minute on 32 CPUs**, without V/J assignment. These differently provisioned measurements do not establish a universal wall-clock or cost winner.[^anarcii-paper]

Its difficult “no-truth” dataset deliberately enriches disagreements between ANARCI reference versions; rates from that set should not be extrapolated to ordinary complete therapeutic antibodies. Some CDR2/DE-region assignments remained ambiguous.[^anarcii-paper]

For deployment, the verified PyPI release is **2.0.8**, dated **30 June 2026**, with **Python ≥3.11** specified. The inspected repository license is **BSD 3-Clause**. The project provides a Python package and documented model-based numbering workflow.[^anarcii-release][^anarcii-license][^anarcii-repo]

My recommendation is to use ANARCII as a serious difficult-case comparator or secondary pass, and evaluate it as a primary engine when GPU batch execution or unusual formats are central. Do not turn a score into a probability of correctness without a separate calibration study. For multi-domain constructs, test domain recovery explicitly rather than assuming a numbering call must return every antibody-like domain.

### 4.3 Immunum: the main Rust-native candidate

Immunum implements a semi-global Needleman–Wunsch-style alignment against position-specific consensus scoring matrices. The core is Rust, with Python bindings, a Polars plugin, and JavaScript/WASM support. It recognizes antibody heavy, kappa, and lambda chains plus four TCR chain classes. The repository is MIT-licensed.[^immunum-source]

This is an attractive engineering combination for a service or data pipeline: native Rust processing can be used directly, while Python users can stay in a tabular workflow rather than repeatedly calling a Python function for each sequence. Those are integration advantages, not independently measured accuracy advantages.

**A version distinction matters.** The verified PyPI release is **1.2.0, uploaded 4 August 2026**. A later source commit, **28 August 2026**, adds Chothia, Martin, and AHo support and changes the source version to 1.3.0. The newer README describes five schemes, with alternate antibody schemes derived from internal IMGT numbering.[^immunum-release][^immunum-commit][^immunum-source]

Therefore, do not assume an unpinned `pip install immunum` obtains every capability shown in the newest source README. Select the released package or pinned source revision deliberately, and test each required scheme and chain class.

#### Has Immunum surpassed AntPack?

The public evidence does not justify a blanket conclusion. Two maintainer issues are particularly relevant:

| Evidence | What it says | Implication |
|---|---|---|
| Issue #32, “Fix `antpack` parallelization benchmark” | Maintainers identify problematic scaling in the AntPack parallel comparison and propose fixes | The comparative speed chart is not a settled basis for declaring a universal winner |
| Issue #33, “Fix correctness benchmarks” | The current reference uses positions where ANARCI and AntPack agree; independent structural truth is requested | Agreement on a consensus subset is not independent structural accuracy |

Both issues were open in the retrieved records.[^immunum-speed-issue][^immunum-truth-issue] These are transparent limitations, not evidence that Immunum is unusable. My recommendation is **first-choice Rust candidate to benchmark**, not **proven most accurate and fastest antibody numberer**.

### 4.4 RIOT: useful when the task is larger than numbering

RIOT combines nucleotide and amino-acid immunoglobulin annotation with an open germline database and includes established numbering conventions. Its peer-reviewed paper appeared online in late 2024 and in the 2025 *Briefings in Bioinformatics* issue.[^riot]

It is therefore relevant when the required output includes germline annotation or repertoire analysis, not just a numbered variable-domain amino-acid sequence. The paper states noncommercial-use conditions for noncommercial organizations; an open database does not imply that every software component is unrestricted for commercial deployment.[^riot]

I would compare RIOT with a complete annotation pipeline, not rank its total runtime against a pure amino-acid numbering kernel without separating the extra work.

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
| **Immunum 1.2.0** | PyPI release 4 August 2026 | Do not infer newer source-only scheme support from the package name |
| **Immunum source `e027d5f…`** | 28 August 2026 commit adds Chothia/Martin/AHo and updates version to 1.3.0 | A pinned source build may be needed for those capabilities at this cutoff |

These records were checked against PyPI and repository sources, not inferred from version ordering alone.[^antpack-current][^antpack-yanked][^antpack-gpl][^anarcii-release][^anarcii-license][^immunum-release][^immunum-commit]

Licensing here is a report of the maintainers' published terms, not legal advice about a specific deployment. Review the actual package, weights, bundled references, and dependencies. In particular, **do not describe all AntPack versions as GPL or all current AntPack versions as commercially unrestricted**.

## 6. Python and Rust examples

These are interface examples grounded in the inspected documentation/source. The engines were not installed or executed in this research session.

### 6.1 Immunum in Python

```python
from immunum import Annotator

# IMGT is available in the verified released interface.
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

For tabular processing, the documented Polars plugin exposes `imp.number(...)` and `imp.segment(...)` as expressions. The same scheme/version warning applies to both the ordinary Python API and the plugin.[^immunum-source][^immunum-release]

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

Use the crate version or source revision corresponding to the API under test. The constructor and method pattern follow the inspected README; the example was not compiled here. Do not assume a stale dependency-version snippet in a README describes the latest released Python and Rust artifacts identically.[^immunum-source]

### 6.3 ANARCII in Python

```python
from anarcii import Anarcii

# Reuse one model instance for a batch rather than constructing one per chain.
model = Anarcii()
imgt_results = model.number([sequence])
martin_results = model.to_scheme("martin")
```

This numbering/conversion pattern is also used by the inspected AbNumber backend. Select and record model mode, device, batching, and version according to the deployed ANARCII documentation; the example intentionally does not imply a particular throughput.[^abnumber][^anarcii-repo]

## 7. Recommended architecture for a semi-automated antibody workflow

The following is a proposed design, not a measured guarantee that a two-engine system is more accurate than either engine alone.

### 7.1 Store a lossless coordinate crosswalk

For each input residue, retain the input sequence offset, domain identity, chain type, scheme, position number, insertion code, and any corresponding observed-structure residue identifier. Store the original sequence and the numbered domain's start/end offsets separately. The displayed label is not an adequate database key on its own.

When processing a structure, number the intended complete domain sequence and map observed coordinates onto it. Missing crystallographic residues should remain explicitly missing; deleting them before numbering can create an artificial shortened sequence. Keep chain identifiers, author residue numbers, insertion codes, and internal sequence offsets distinct.

This crosswalk supports comparisons between schemes without destroying the original coordinates. It also lets the superposition workflow select a framework for fitting and a CDR for evaluation without confusing either with arbitrary MSA column numbers.

### 7.2 Use a fast primary pass and explicit escalation criteria

For a permissively licensed Rust/Python service, I would begin by evaluating Immunum against an internal truth panel. For a CPU-first Python research pipeline, I would evaluate a license-compatible AntPack version as the throughput baseline. Use ANARCII as an independent difficult-case comparator or as the primary engine for workloads where its strengths are important.

Escalation criteria should include domain-detection failure, unusually low score, missing or duplicated labels, unexpected truncation, unusual insertions, disagreement in CDR boundaries, long multi-domain constructs, and structurally inconsistent anchor placement. Agreement between two engines is reassuring only to the extent that their errors are independent; do not equate consensus with truth.

### 7.3 Keep numbering and biological annotations modular

Numbering, CDR segmentation, germline assignment, humanness, and developability assessment should be distinct versioned stages. A tool that performs several of these may be convenient, but compare its numbering performance with optional analyses disabled before making a kernel-speed claim.

Convert schemes only after preserving the original mapping. An IMGT-to-Martin conversion is not evidence that the input was newly aligned against a Martin-specific structural reference. That distinction is especially important when inspecting disagreements between engines that share an IMGT intermediate.[^immunum-source][^abnumber]

## 8. Benchmark that would answer “best for us”

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

## 9. Bottom line

**Scheme choice:** Support IMGT for interchange and retain Martin or AHo when structural or legacy analyses require them. A new numbering engine is not, by itself, a new and structurally superior scheme.

**Engine choice:** AntPack remains a strong CPU baseline; ANARCII has meaningful recent evidence on difficult sequence classes; Immunum is the most directly relevant Rust-native candidate found in this research. RIOT and AbNumber solve useful adjacent workflow problems, but their scope and backend should be made explicit.

**Has AntPack been surpassed?** On particular difficult accuracy tests, there is published evidence in favor of ANARCII. For a universal speed-and-accuracy claim, or a settled claim that Immunum is faster on a fair matched benchmark, the retrieved evidence is insufficient. The appropriate next decision is a version-pinned benchmark on the project's sequence distribution, not an unconditional leaderboard.

## References and implementation records

[^anarci]: Dunbar J, Deane CM. **ANARCI: antigen receptor numbering and receptor classification.** *Bioinformatics* 32, 298–300 (2016; online 2015). DOI: <https://doi.org/10.1093/bioinformatics/btv552>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/26424857/>.
[^immunum-source]: ENPICOM, **Immunum README**, including algorithm, interfaces, chain/scheme support, and MIT license declaration: <https://github.com/ENPICOM/immunum/blob/e027d5fe2405300508eee7ce78582a2fe4990f20/README.md>.
[^martin]: Abhinandan KR, Martin ACR. **Analysis and improvements to Kabat and structurally correct numbering of antibody variable domains.** *Molecular Immunology* 45, 3832–3839 (2008). DOI: <https://doi.org/10.1016/j.molimm.2008.05.022>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/18614234/>.
[^aho]: Honegger A, Plückthun A. **Yet another numbering scheme for immunoglobulin variable domains: an automatic modeling and analysis tool.** *Journal of Molecular Biology* 309, 657–670 (2001). DOI: <https://doi.org/10.1006/jmbi.2001.4662>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/11397087/>.
[^evaluation2024]: Zhu Z, Olson KS, Magliery TJ. **50 Years of Antibody Numbering Schemes: A Statistical and Structural Evaluation Reveals Key Differences and Limitations.** *Antibodies* 13, 99 (2024). DOI: <https://doi.org/10.3390/antib13040099>. Publisher: <https://www.mdpi.com/2073-4468/13/4/99>; PubMed: <https://pubmed.ncbi.nlm.nih.gov/39727482/>.
[^antpack-paper]: Parkinson J, Wang W. **For antibody sequence generative modeling, mixture models may be all you need.** *Bioinformatics* 40, btae278 (2024). DOI: <https://doi.org/10.1093/bioinformatics/btae278>. This paper includes the AntPack numbering benchmark: <https://academic.oup.com/bioinformatics/article/40/5/btae278/7656770>.
[^anarcii-paper]: Greenshields-Watson A, Agarwal P, Robinson SA et al. **ANARCII enables alignment-free antigen receptor numbering using a generalised language model.** *Communications Biology* 9, 1085 (2026). Published 21 May 2026; version-of-record date 12 August 2026. DOI: <https://doi.org/10.1038/s42003-026-10186-z>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/42162238/>. Numerical comparisons in this report are author-reported, not reproduced here.
[^riot]: Dudzic et al. **RIOT—Rapid Immunoglobulin Overview Tool—annotation of nucleotide and amino acid immunoglobulin sequences using an open germline database.** *Briefings in Bioinformatics* 26, bbae632 (2025 issue; online 2024). DOI: <https://doi.org/10.1093/bib/bbae632>. Publisher: <https://academic.oup.com/bib/article/26/1/bbae632/7914577>; official project: <https://github.com/NaturalAntibody/riot_na>.
[^immunum-speed-issue]: Immunum maintainer issue **#32**, *Fix `antpack` parallelization benchmark*: <https://github.com/ENPICOM/immunum/issues/32>. Opened 23 March 2026, updated 8 April 2026; open in the retrieved record.
[^immunum-truth-issue]: Immunum maintainer issue **#33**, *Fix correctness benchmarks*: <https://github.com/ENPICOM/immunum/issues/33>. Opened 23 March 2026; open in the retrieved record. The issue explicitly requests an independent, for example structure-based, gold standard.
[^antpack-current]: AntPack, current PyPI package metadata and licensing: <https://pypi.org/project/antpack/>. At retrieval, the default non-yanked release was 0.4; versions 0.3.9 onward use academic/noncommercial terms and key setup.
[^antpack-gpl]: AntPack **0.3.8.6.3**, release record and GPL metadata: <https://pypi.org/project/antpack/0.3.8.6.3/>. Released 23 June 2026; Python ≥3.8 stated in this package record.
[^immunum-release]: Immunum **1.2.0**, PyPI package and release history: <https://pypi.org/project/immunum/>. Source distribution and wheels uploaded 4 August 2026.
[^immunum-commit]: Immunum commit **`e027d5fe2405300508eee7ce78582a2fe4990f20`**, 28 August 2026, adding Chothia, Martin, and AHo and changing source version to 1.3.0: <https://github.com/ENPICOM/immunum/commit/e027d5fe2405300508eee7ce78582a2fe4990f20>.
[^antpack-regions]: AntPack, official **clustering / region assignment** documentation, including separate numbering and CDR conventions: <https://antpackdocumentationlatest.pages.dev/clustering_overview>.
[^imgt]: IMGT, **IMGT unique numbering for V domains**, official scientific chart: <https://imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html>.
[^antpack-doc]: AntPack, official **numbering background** documentation: <https://antpackdocumentationlatest.pages.dev/numbering_background>. The documentation site's displayed version can lag package releases; verify the deployed API separately.
[^anarcii-repo]: ANARCII, inspected release README and official user guide: <https://github.com/oxpig/ANARCII/blob/e0d8f192f5a861e03a50918f114d0f5735e42333/README.md>; <https://github.com/oxpig/ANARCII/wiki>.
[^abnumber]: AbNumber, inspected **`abnumber/common.py`**, `_anarci_align` backend selection, scheme conversion, germline fallback, and duplicate-position handling: <https://github.com/prihoda/AbNumber/blob/master/abnumber/common.py>. File blob hash at retrieval: `d494cb3593745ee824727537176da5eec67e1aaf`. This is an inspected source record, not a guarantee about every packaged version.
[^anarcii-release]: ANARCII **2.0.8** package metadata and release history: <https://pypi.org/project/anarcii/>. Release date 30 June 2026; Python ≥3.11 in the retrieved package metadata.
[^anarcii-license]: ANARCII, pinned **BSD 3-Clause** license: <https://github.com/oxpig/ANARCII/blob/e0d8f192f5a861e03a50918f114d0f5735e42333/LICENCE>.
[^antpack-yanked]: AntPack **0.5** release record and release history: <https://pypi.org/project/antpack/0.5/>. Released 14 April 2026; yanked for “Bug fix.”
