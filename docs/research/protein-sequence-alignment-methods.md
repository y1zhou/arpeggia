# Pairwise and multiple protein sequence alignment: accepted methods and developments through September 2026

**Research cutoff:** 8 September 2026  
**Scope:** Amino-acid sequence alignment, with distinctions between exact pairwise optimization, database search, multiple alignment, and structure-informed correspondence.  
**Evidence:** Primary publications and official software documentation. Recommendations below are a task-based synthesis, not a measured ranking of software market share.  
**Execution:** The core Biopython example was checked with Biopython 1.86. No comparative performance benchmark was run in this session.

## Executive answer

**Smith–Waterman has not become obsolete.** For a specified local-alignment scoring model, an exact implementation already finds an optimal score. A newer algorithm cannot obtain a strictly better optimum for that identical problem. It can run faster, consume less memory, return a different equally optimal path, or use a richer model that better reflects biological correspondence.[^sw][^gotoh]

For two known protein sequences, my starting choices are **affine-gap global or local dynamic programming**, exposed through Biopython for convenience or a SIMD implementation such as Parasail for larger workloads. Choose the alignment objective before selecting the fastest implementation.[^biopython][^parasail]

For protein MSAs, **MAFFT, MUSCLE 5, and Clustal Omega** remain useful established reference choices, but the appropriate mode and dataset size matter. **FAMSA2**, published in April 2026, belongs on the current shortlist for large protein families; **UPP2** is particularly relevant to fragment-rich datasets.[^mafft][^muscle5][^clustalo][^famsa2][^upp2]

The important recent advances are not one universal replacement: they include faster exact kernels, more selective database-search pipelines, better scalable MSA construction, and learned or structural information for difficult correspondences.

## 1. First identify the problem being solved

| Task | Appropriate starting point | What would be a category mistake? |
|---|---|---|
| Two homologous, full-length domains with comparable architecture | Global alignment with a protein substitution matrix and affine gaps | Using a local score that silently ignores poorly matching ends |
| One shared domain inside otherwise different proteins | Local alignment, or detect and align the domains separately | Forcing unrelated domains into one end-to-end alignment |
| A fragment versus a full domain | Explicit overlap / semi-global end-gap policy | Penalizing missing termini as if they were internal evolutionary deletions |
| One sequence versus millions of database entries | BLAST/MMseqs2-style search pipeline; profile search when appropriate | Running a slow Python pairwise loop over the whole database |
| A family of homologous proteins | An MSA method appropriate to family size and length variation | Assuming pairwise optimality implies an optimal MSA |
| Remote relationships with structures available | Structure-informed cross-check or alignment | Calling improved structural correspondence a faster solution to the same sequence-only objective |
| Antibody positional annotation | A numbering engine and explicit scheme | Treating arbitrary MSA columns as interchangeable with IMGT or AHo positions |

The first three distinctions correspond to the local/global/end-gap options documented for pairwise alignment; the database-search and MSA rows describe different algorithmic pipelines.[^biopython][^mmseqsgpu][^mafft][^upp2] The antibody distinction is examined in the companion numbering report.

## 2. Classical pairwise alignment is still the foundation

### 2.1 Local, global, and semi-global alignment

**Global alignment**, associated with Needleman–Wunsch, aligns the sequences end to end. **Local alignment**, associated with Smith–Waterman, identifies the best-scoring subsequences. **Semi-global or overlap alignment** changes which terminal gaps are penalized. “Semi-global” is not a complete specification by itself: all four sequence ends can have distinct policies.[^biopython][^sw]

For protein design, this choice affects interpretation. A high-scoring local match may cover only an intact framework while excluding a redesigned loop or terminal extension. That can be useful for detecting a domain but inappropriate for claiming that the complete design closely matches the reference.

### 2.2 Affine gaps and substitution matrices

A common scoring objective is

$$
S=\sum_{(i,j)\in A}M(a_i,b_j)
-\sum_{g\in G}\left[g_o+(|g|-1)g_e\right],
$$

where $M$ is a protein substitution matrix, $A$ contains aligned residue pairs, and each gap run has an opening and extension penalty. Gotoh's algorithm organizes affine-gap dynamic programming efficiently using distinct gap states.[^gotoh]

BLOSUM62 is a reasonable starting matrix for many protein comparisons, not a universal biological truth. Divergence, domain architecture, repeats, low-complexity regions, and nonstandard residues can all make a default alignment misleading. The correct output is a model-dependent correspondence, not proof of homology merely because an algorithm can print an alignment.[^biopython][^dedal]

A reproducibility trap is the gap convention: some APIs use $g_o+(L-1)g_e$, others use $g_o+Lg_e$. Transferring identical numeric parameters between these conventions changes the score of every gap run. Also distinguish libraries accepting negative gap scores from those accepting positive costs.[^biopython][^wfa]

### 2.3 What “optimal” does and does not mean

For lengths $m,n$, conventional exact affine-gap DP requires quadratic $O(mn)$ time. Score-only calculation can keep a small number of rows; full traceback generally requires more storage unless an explicit low-memory traceback strategy is used. Returning only a score and returning a full residue mapping are different workloads.[^gotoh][^biopython][^biwfa]

Exact score optimization cannot recover information that is absent from its scoring model. Two low-identity sequences may have several biologically plausible alignments, and repeats can produce equally scoring alternatives. A richer model can improve correspondence accuracy while still using dynamic programming as the optimizer. That is a change in the inference model, not a refutation of Smith–Waterman's optimality.[^dedal]

## 3. What has improved runtime?

### 3.1 Faster exact implementations: SIMD and careful engineering

Parasail implements vectorized local, global, and semi-global pairwise alignment. Its SIMD kernels calculate multiple DP operations in parallel without inherently changing the chosen scoring objective. It is a useful option when many already-selected protein pairs must be aligned and Python convenience remains desirable.[^parasail]

My practical advice is to first remove avoidable overhead: construct reusable scoring objects once, batch work, avoid subprocess startup for every short sequence, and use score-only kernels when traceback is unnecessary. For compact integer SIMD implementations, check score saturation and the available wider-integer fallback. A saturated score is not an accurate result merely because the computation completed.[^parasail]

### 3.2 Heuristic acceleration for protein scores: Block Aligner

Block Aligner is a Rust implementation designed around adaptively shifted and resized blocks of the alignment matrix. It supports protein-relevant scoring and C bindings. Unlike an unrestricted exact DP kernel, the adaptive search is heuristic; the paper reports substantial speed gains with a small but nonzero error rate on its tested workloads.[^block]

This makes it a strong candidate for a Rust-heavy, high-volume pipeline where occasional approximate results are acceptable or can be rechecked. It is not my first choice for a small, correctness-critical mutation-mapping task where exact DP is already inexpensive. That recommendation follows from the different failure costs, rather than a claim that the method is universally inferior or superior.

### 3.3 Wavefront alignment: an important advance with a protein-specific caveat

The wavefront algorithm, WFA, gives exact gap-affine alignment for its supported edit-style scoring model. Its cost depends on the optimal alignment score $s$, yielding $O(ns)$ time and $O(s^2)$ working space in the original formulation. It is especially attractive when sequences are similar enough that $s$ is small.[^wfa]

BiWFA reduces working memory to $O(s)$ while retaining the favorable score-dependent time bound. These are meaningful algorithmic advances, not merely faster code.[^biwfa]

**However, the standard WFA formulation is not a drop-in replacement for arbitrary BLOSUM/PAM-scored protein Smith–Waterman.** Uniform mismatch/edit-style costs do not express a general amino-acid substitution matrix. Worst-case behavior also loses the attractive near-linear appearance when divergence makes $s$ large.[^wfa][^biwfa]

A useful decision rule is therefore: use WFA-family methods when the actual scoring and end-gap objective match the implementation; do not select them solely because the input alphabet happens to contain amino-acid letters. Similarly, a very fast unit-cost edit-distance library does not automatically provide an appropriate biological protein alignment.

### 3.4 Database search: avoid most full alignments

Database-search tools obtain much of their speed by identifying promising candidates before doing expensive gapped alignment. This changes the end-to-end search strategy, not necessarily the final alignment kernel. Search sensitivity therefore needs its own benchmark; exact alignment of shortlisted candidates does not prove that every relevant candidate was shortlisted.[^mmseqsgpu]

The 2025 MMseqs2-GPU paper is a useful contemporary example. It combines GPU-accelerated filtering with a gapped Smith–Waterman–Gotoh stage. Thus a modern, much faster search system can still rely on the classical exact recurrence in its final alignment component.[^mmseqsgpu]

For structure-prediction pipelines, “build an MSA” often means search, filter, cluster, and assemble homologs around a query. This is not identical to aligning a fixed, already-collected family with MAFFT or FAMSA2. Benchmark the complete intended pipeline rather than comparing these different tasks as though they were interchangeable.

## 4. What has improved biological alignment accuracy?

### 4.1 More information than a fixed pair of sequences

Profile-based methods represent position-specific residue preferences and gap behavior. An MSA can therefore make a difficult pairwise correspondence easier by contributing information from other homologs. This also explains why profile-oriented MSA methods can outperform repeated independent pairwise alignment on a divergent family.[^clustalo][^upp2]

The cost is that profile construction and family selection become part of the inference. Incorrectly grouped domains or overrepresented sequence subfamilies can bias the result. My recommendation is to assess domain composition and sequence diversity before interpreting a profile-derived alignment as a definitive residue map.

### 4.2 Learned contextual scoring

DEDAL, published online in December 2022 and in the 2023 issue of *Nature Methods*, learns protein representations and alignment-relevant scores rather than relying solely on a fixed substitution matrix. Its authors reported marked improvements for remote homologs. This is a concrete example of improved biological correspondence through richer scoring information.[^dedal]

The engineering question is whether the accuracy gain on the intended sequence class justifies model inference, hardware, and deployment costs. A language-model-derived alignment is not automatically faster than classical DP, and its evaluation must control for training overlap and difficult out-of-distribution examples.

### 4.3 Structural evidence

When reliable structures are available, structural alignment provides a separate source of positional evidence. FoldMason's peer-reviewed 2026 work addresses scalable multiple protein structure alignment. It is relevant as a structural comparison or validation layer, not as a sequence-only method that optimizes exactly the same objective as Smith–Waterman.[^foldmason]

For designed proteins, structural and evolutionary correspondence can disagree for legitimate reasons. A redesigned loop may preserve geometry without preserving sequence homology, while a homologous flexible region can adopt a different conformation. State which notion of equivalence the analysis is intended to recover.

## 5. Multiple sequence alignment: the practical landscape

MSA programs usually combine a guide tree, progressive sequence/profile alignment, and sometimes consistency transformations or iterative refinement. They do not generally guarantee the globally best alignment under every possible multi-sequence objective. Early alignment decisions, guide-tree quality, and how fragments are handled can affect the final result.[^tcoffee][^mafft][^muscle5]

### 5.1 Current decision table

| Method | Main role in a 2026 workflow | Accuracy / runtime considerations |
|---|---|---|
| **MAFFT L-INS-i / G-INS-i / E-INS-i** | High-quality small-to-medium protein-family alignment with selectable assumptions | Match the mode to local homology, full-length homology, or large unalignable inserts |
| **MAFFT faster / automatic modes** | Practical general-purpose starting point | `--auto` is a strategy selector, not a certificate of optimal accuracy |
| **MUSCLE 5** | Accuracy-oriented alignment and alignment ensembles | Useful when downstream conclusions depend on alignment uncertainty; do not confuse it with old MUSCLE versions |
| **Clustal Omega** | Established scalable profile-based baseline | Useful for interoperability and comparative benchmarking, not necessarily the fastest current option |
| **T-Coffee family** | Consistency-based alignment / combination of alignment evidence | Established accuracy-oriented approach; cost can be higher than simpler progressive methods |
| **FAMSA2** | Large protein families and high-throughput MSA | Important 2026 addition to the shortlist; use the actual new backend |
| **UPP2** | Many fragments and strong sequence-length heterogeneity | Backbone alignment plus profile-HMM placement addresses a different problem from uniform full-length families |
| **ARIES** | Emerging embedding-based MSA, particularly difficult low-identity families | Promising 2026 research; include embedding-generation cost and test independently |
| **FoldMason** | Structure-informed multiple alignment | Requires structural information and should be evaluated as a different input regime |

The method descriptions are grounded in the primary papers and official MAFFT documentation.[^mafft][^mafft-doc][^muscle5][^clustalo][^tcoffee][^famsa2][^upp2][^aries][^foldmason]

### 5.2 MAFFT: choose the actual mode

The relevant distinctions are L-INS-i for local pairwise information and iterative refinement, G-INS-i for globally alignable sequences, and E-INS-i for sequences containing large unalignable regions. Fast MAFFT modes make different accuracy–cost tradeoffs. For fragments, the documented addition modes can be more appropriate than rebuilding an MSA while treating every sequence as full length.[^mafft][^mafft-doc]

For an ordinary collection of complete homologous domains, I would use a suitable accuracy-oriented MAFFT mode as one reference alignment. For large collections, I would compare it on a representative subset and use a scalable full-dataset method rather than assuming a costly mode must remain best at every scale.

### 5.3 MUSCLE 5: alignment uncertainty is a first-class output

MUSCLE 5 emphasizes high-accuracy alignment and ensembles generated through perturbations of alignment construction. The point is not only obtaining one good MSA; it is assessing how conclusions change across plausible alternatives. Its scalable Super5 mode is another reason to specify the exact version and mode, not merely “MUSCLE.”[^muscle5]

For sensitive evolutionary or positional analyses, compare downstream conclusions across an alignment ensemble or across well-motivated methods. Bootstrap resampling of one fixed MSA does not by itself explore uncertainty in the MSA's residue correspondences. This is an analysis recommendation motivated by the ensemble approach.

### 5.4 FAMSA2: a material update in April 2026

FAMSA2 was published in *Nature Biotechnology* on **14 April 2026**. It combines progressive alignment with medoid-based guide-tree construction and an LCS-derived dissimilarity measure. The authors report matching or exceeding competing accuracy across structural, phylogenetic, and functional benchmarks, with an average **400-fold runtime advantage** in their comparisons.[^famsa2]

That headline is a study-level result, not a prediction of a 400-fold gain on every family or machine. The accessible publication page exposed the abstract, methods overview, and extended-data descriptions; the complete subscription main text was not independently audited here. My recommendation is to benchmark FAMSA2 for large families and retain an accuracy-focused comparator on a curated subset. Also verify that a Python wrapper actually embeds FAMSA2 rather than an older FAMSA release.

### 5.5 UPP2: length heterogeneity deserves its own solution

UPP2 aligns a backbone and uses an ensemble of profile HMMs to place remaining sequences, improving the efficiency of the earlier UPP strategy. Its focus on fragmentary data makes it especially relevant when length distributions are broad.[^upp2]

For such datasets, a benchmark containing only complete proteins can select the wrong tool. Measure whether fragments land in the correct domain and whether they distort the full-length backbone; one aggregate column score can conceal both problems.

### 5.6 ARIES: promising, but not yet a universal production default

The 2026 ARIES work uses protein-language embeddings, a reciprocal embedding similarity measure, and a template-based dynamic-time-warping construction. Its authors report improved low-identity accuracy and favorable scaling. The accessible primary record was a 2026 preprint, so it is treated here as emerging evidence rather than an established replacement for all classical MSA tools.[^aries]

My recommendation is to include it in a difficult-family pilot, with embeddings charged to the runtime budget. Cached-embedding timing and end-to-end timing answer different questions. Its value should be tested on the intended family distribution, not inferred from publication recency alone.

## 6. Python and Rust integration

### 6.1 Python: an explicit, auditable pairwise baseline

Biopython's `Bio.Align.PairwiseAligner` is the current documented interface to use for this purpose. It exposes substitution matrices, gap scoring, alignment modes, and coordinate mappings.[^biopython]

```python
from Bio.Align import PairwiseAligner, substitution_matrices


def align_proteins(seq_a: str, seq_b: str, mode: str = "global"):
    """Return one optimal alignment and its zero-based ungapped residue pairs.

    Input is an unaligned amino-acid sequence; stop symbols and gap characters
    are deliberately rejected. Other alphabets require an explicit policy.
    """
    if mode not in {"global", "local"}:
        raise ValueError("mode must be 'global' or 'local'")

    matrix = substitution_matrices.load("BLOSUM62")
    allowed = set(matrix.alphabet) - {"*"}
    seq_a, seq_b = seq_a.upper().strip(), seq_b.upper().strip()
    for label, sequence in (("seq_a", seq_a), ("seq_b", seq_b)):
        if not sequence:
            raise ValueError(f"{label} is empty")
        invalid = set(sequence) - allowed
        if invalid:
            raise ValueError(f"Unsupported residues in {label}: {sorted(invalid)}")

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq_a, seq_b)
    try:
        alignment = alignments[0]
    except IndexError as exc:
        raise ValueError("No positive-scoring local alignment was found") from exc

    # -1 indicates a gap; this mapping is independent of display formatting.
    pairs = [(int(i), int(j)) for i, j in alignment.indices.T
             if i >= 0 and j >= 0]
    return alignment, pairs


alignment, residue_pairs = align_proteins("ACDEFGHIK", "ACDEYGHIK")
print(alignment.score)  # 50.0 with these explicit parameters
print(residue_pairs)
```

The underlying construction, score, and `.indices` interface were exercised in Biopython 1.86. The penalties are an illustrative explicit baseline, not an asserted biological optimum. For high-volume work, move matrix and aligner construction outside the per-sequence call. When only a score is required, use the score-only interface instead of constructing the traceback.[^biopython]

Do not enumerate every equally optimal alignment without considering the possible output size. For a mutation-coordinate pipeline, store the chosen mapping and tie policy rather than relying on an alignment's printed string.

### 6.2 Python throughput and Rust-native use

| Requirement | Candidate | Qualification |
|---|---|---|
| Simple Python pairwise mapping | Biopython `PairwiseAligner` | Explicit mode and gap settings; convenient reference implementation |
| Many exact protein pairs | Parasail with Python bindings | Verify traceback, integer width, saturation, and end-gap mode |
| Rust-native fast approximate protein alignment | Block Aligner | Benchmark heuristic error and fall back to exact DP where needed |
| Rust integration with an exact native kernel | A tested C ABI to a mature alignment library | Extra build/FFI work; preserve score and traceback conventions |
| Protein MSA from Python | Invoke a pinned MAFFT, MUSCLE 5, FAMSA2, or UPP2 backend | Wrapper convenience does not determine the underlying algorithm/version |

These options follow the interfaces and algorithms described in their primary documentation and papers.[^biopython][^parasail][^block][^mafft-doc][^muscle5][^famsa2][^upp2] For a Rust-based antibody pipeline specifically, Immunum is discussed in the separate numbering report; it is not a general replacement for arbitrary protein pairwise alignment.

Example MAFFT invocations from its documented modes are:

```bash
# Complete domains with local conserved blocks: L-INS-i.
mafft --localpair --maxiterate 1000 proteins.fasta > proteins.linsi.fasta

# Globally alignable complete domains: G-INS-i.
mafft --globalpair --maxiterate 1000 proteins.fasta > proteins.ginsi.fasta

# Large unalignable insertions: E-INS-i.
mafft --genafpair --maxiterate 1000 proteins.fasta > proteins.einsi.fasta
```

These commands were not executed in this session. Pin executable versions and retain stderr logs, parameters, and input ordering.[^mafft-doc]

## 7. How to compare tools without misleading yourself

The following is a proposed evaluation protocol, not a benchmark performed for this report.

### 7.1 Pairwise validation

For an exact-kernel comparison, fix the substitution matrix, gap convention, local/global/end-gap objective, traceback requirements, and treatment of ambiguous residues. Compare optimal **scores** first. Different paths with equal scores are not necessarily implementation errors; residue-map differences still matter to downstream applications.

Then test representative failure cases: complete domains, fragments, long insertions, repeats, low identity, low complexity, and designed sequences. Report sequence coverage alongside identity and score. Define the identity denominator explicitly—aligned residue pairs, alignment columns, shorter sequence, and query length are not interchangeable.

Separate algorithm time from parsing, object construction, thread startup, process startup, data transfer, and output serialization. Report both cold and warm runs. Score-only CPU timing must not be compared unqualified with full-traceback GPU timing.

### 7.2 MSA validation

Use at least one benchmark with independently curated or structure-supported correspondences, and separate family classes rather than relying on a single overall average. Sum-of-pairs recovery and complete-column recovery measure different aspects of accuracy; report coverage and runtime/memory failures as well.

For a protein-design application, add application-level tests: consistency of conserved motif placement, stability of mutation-coordinate maps, conservation estimates at selected sites, and agreement with trustworthy structural anchors. An alignment can score well on average while misplacing the one loop or interface segment that matters to a design decision.

For learned methods, audit train/test overlap and include difficult held-out families. For predicted-structure evaluation, identify where structure confidence is low and avoid treating uncertain coordinates as unquestionable ground truth. For multi-domain inputs, evaluate domain matching separately from within-domain alignment.

## 8. Recommended stack for this project

My proposed default architecture is an **explicit correspondence layer**, not one command used for every problem.

For two known related domains, start with an exact affine-gap mapping. Use global or overlap alignment when the intended correspondence spans the domain; use local alignment for domain discovery or partial homology. Keep a record of the score model and the selected residue pairs.

For family-level analyses, use an accuracy-oriented MAFFT or MUSCLE 5 run as a reference on manageable subsets, assess FAMSA2 for large-scale production, and introduce UPP2 when fragments dominate. Examine alternative alignments in regions where a design conclusion depends on an ambiguous correspondence.

For remote or structurally unusual sequences, add contextual or structure-informed evidence rather than trying to solve missing biological information solely by tuning gap penalties. For antibodies, introduce a numbering-aware mapping before deriving CDR mutation coordinates or structural evaluation selections.

The overarching answer is therefore **yes, there have been substantial improvements in runtime and biological alignment quality—but no single algorithm “surpasses Smith–Waterman” across all these different objectives and input regimes.**

## References

[^sw]: Smith TF, Waterman MS. **Identification of common molecular subsequences.** *Journal of Molecular Biology* 147, 195–197 (1981). DOI: <https://doi.org/10.1016/0022-2836(81)90087-5>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/7265238/>.
[^gotoh]: Gotoh O. **An improved algorithm for matching biological sequences.** *Journal of Molecular Biology* 162, 705–708 (1982). DOI: <https://doi.org/10.1016/0022-2836(82)90398-9>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/7166760/>.
[^biopython]: Biopython, **Pairwise sequence alignment**, official tutorial: <https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html>. Documentation checked at the research cutoff; local API example checked with Biopython 1.86.
[^parasail]: Daily J. **Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments.** *BMC Bioinformatics* 17, 81 (2016). DOI: <https://doi.org/10.1186/s12859-016-0930-z>. Official implementation: <https://github.com/jeffdaily/parasail>.
[^block]: Liu D, Steinegger M. **Block Aligner: an adaptive SIMD-accelerated aligner for sequences and position-specific scoring matrices.** *Bioinformatics* 39, btad487 (2023). DOI: <https://doi.org/10.1093/bioinformatics/btad487>. Full text: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10457662/>.
[^wfa]: Marco-Sola S et al. **Fast gap-affine pairwise alignment using the wavefront algorithm.** *Bioinformatics* (2021). Full text: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8355039/>.
[^biwfa]: Marco-Sola S et al. **Optimal gap-affine alignment in O(s) space.** *Bioinformatics* 39, btad074 (2023). DOI: <https://doi.org/10.1093/bioinformatics/btad074>. Full text: <https://pmc.ncbi.nlm.nih.gov/articles/PMC9940620/>.
[^mmseqsgpu]: **GPU-accelerated homology search with MMseqs2.** *Nature Methods* 22, 2024–2027 (2025). DOI: <https://doi.org/10.1038/s41592-025-02819-8>.
[^dedal]: Llinares-López F et al. **Deep embedding and alignment of protein sequences.** *Nature Methods* 20, 104–111 (2023; online 15 December 2022). DOI: <https://doi.org/10.1038/s41592-022-01700-2>.
[^mafft]: Katoh K, Standley DM. **MAFFT multiple sequence alignment software version 7: improvements in performance and usability.** *Molecular Biology and Evolution* 30, 772–780 (2013). DOI: <https://doi.org/10.1093/molbev/mst010>.
[^mafft-doc]: MAFFT, **official software site and usage documentation**: <https://mafft.cbrc.jp/alignment/software/>; <https://mafft.cbrc.jp/alignment/software/manual/manual.html>.
[^muscle5]: Edgar RC. **Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny.** *Nature Communications* 13, 6968 (2022). DOI: <https://doi.org/10.1038/s41467-022-34630-w>.
[^clustalo]: Sievers F et al. **Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega.** *Molecular Systems Biology* 7, 539 (2011). DOI: <https://doi.org/10.1038/msb.2011.75>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/21988835/>.
[^tcoffee]: Notredame C, Higgins DG, Heringa J. **T-Coffee: A novel method for fast and accurate multiple sequence alignment.** *Journal of Molecular Biology* 302, 205–217 (2000). DOI: <https://doi.org/10.1006/jmbi.2000.4042>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/10964570/>.
[^famsa2]: Gudyś A, Zielezinski A, Notredame C, Deorowicz S. **Fast and accurate multiple-protein-sequence alignment at scale with FAMSA2.** *Nature Biotechnology*, published 14 April 2026. DOI: <https://doi.org/10.1038/s41587-026-03095-3>. Official implementation: <https://github.com/refresh-bio/FAMSA>. The quoted speedup is author-reported and dataset-dependent.
[^upp2]: **UPP2: fast and accurate alignment of datasets with fragmentary sequences.** *Bioinformatics* 39, btad007 (2023). DOI: <https://doi.org/10.1093/bioinformatics/btad007>.
[^aries]: Hoang M, Armour-Garb I, Singh M. **Fast, accurate construction of multiple sequence alignments from protein language embeddings.** ARIES; 2026 preprint. DOI: <https://doi.org/10.64898/2026.01.02.697423>. Primary record: <https://www.biorxiv.org/content/10.64898/2026.01.02.697423v1>; full-text record: <https://pmc.ncbi.nlm.nih.gov/articles/PMC13060855/>. First posted 2 January 2026; the indexed full-text record includes a later March revision.
[^foldmason]: Gilchrist CLM, Mirdita M, Steinegger M. **Multiple protein structure alignment at scale with FoldMason.** *Science* 391, 485–488 (2026). DOI: <https://doi.org/10.1126/science.ads6733>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/41610233/>.
