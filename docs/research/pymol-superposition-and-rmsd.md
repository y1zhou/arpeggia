# PyMOL superposition, iterative atom rejection, and RMSD

**Research cutoff:** 8 September 2026
**Scope:** Open-source PyMOL; source audit plus scientific context.
**Audited revision:** `5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69` (24 July 2026), the latest default-branch commit returned by the repository API before the cutoff.[^revision]
**Validation status:** The relevant Python and C++ source paths were inspected. PyMOL itself was not installed or executed in this research session. An independent NumPy demonstration of the rejection rule was executed; it is not a binary-level PyMOL regression test.

## Executive findings

The most important finding is a documentation–implementation discrepancy: **in the audited implementation, `align(..., cutoff=2.0)` rejects an atom pair when its post-fit distance exceeds twice the current RMSD—not when it exceeds 2 Å.** The Python docstring describes an Å cutoff, but the C++ expression divides the distance by RMSD before comparing it with `cutoff`.[^wrapper-align][^rejection]

`align` first constructs a sequence-based residue correspondence and then matching atom pairs. Its refinement loop does **not** repeatedly align the sequences or search for new nearest-neighbor atom correspondences. It repeatedly fits the current pairs, removes outlying pairs, and refits the survivors. Rejected pairs never re-enter that invocation.[^executive][^rejection]

The reported final RMSD is consequently a **surviving-core RMSD**. It is not necessarily the RMSD of all residues, all selected atoms, or all atoms in the two structures. Even `cycles=0` measures only the atom pairs admitted by the initial correspondence.[^selector][^rejection]

For quantitative protein-design evaluation, my recommendation is to retain both an independently defined evaluation correspondence and the fitting correspondence. Report the atom selection, initial and final pair counts, coverage, and RMSD with and without rejection. A low trimmed RMSD by itself is an incomplete comparison.

## 1. Three different operations hidden behind “alignment”

For two coordinate arrays, three decisions should be kept separate:

1. **Correspondence:** Which residue and atom in structure A corresponds to which in B?
2. **Superposition:** Which rigid rotation and translation best fit the chosen corresponding coordinates?
3. **Evaluation:** On which pairs is the final RMSD calculated, and under which transformation?

PyMOL's `align` combines these decisions into a convenient command, but its source implements them in distinct stages.[^executive] The conceptual distinction matters: changing the correspondence is not equivalent to changing the numerical optimizer, and evaluating only the fitted core is not equivalent to evaluating the full protein.

For a set of paired atoms $S$, the equal-weight rigid least-squares objective is

$$
\operatorname{RMSD}(S)=\min_{R,t}\sqrt{\frac{1}{|S|}\sum_{i\in S}\|R x_i+t-y_i\|^2},
\qquad R^TR=I,\quad\det R=1.
$$

Kabsch's work provides the classical mathematical context for this optimization.[^kabsch] However, **“the objective is the Kabsch least-squares problem” does not establish that every PyMOL build executes a conventional SVD implementation of Kabsch**. The audited `MatrixFitRMSTTTf` has alternative numerical branches, and the source default for `fit_kabsch` is zero.[^matrix][^settings]

## 2. The source-level call path

| Stage | Source location / function | What it determines |
|---|---|---|
| Python interface | `modules/pymol/fitting.py`: `align`, `super` | Defaults, selections, states, and arguments forwarded to C++ |
| Command binding | `layer4/Cmd.cpp` | Invocation and assembly of the returned statistics |
| Residue correspondence | `layer3/Executive.cpp`: `ExecutiveAlign`; `layer0/Match.cpp`: `MatchAlign` | Sequence or structural residue-pair scoring and alignment path |
| Atom correspondence | `layer3/Selector.cpp`: `SelectorCreateAlignments` | Selected matching atoms within aligned residue pairs |
| Fit and rejection | `layer3/Executive.cpp`: `ExecutiveRMS` | Initial fit, normalized-distance rejection, repeated fitting |
| Numerical fit | `layer0/Matrix.cpp`: `MatrixFitRMSTTTf` | Rigid-body transformation and RMSD for the current coordinate lists |

These are different layers of one pipeline, not independent external alignment programs.[^wrapper-align][^executive][^match][^selector][^matrix][^binding]

### 2.1 What `align` uses for its initial residue alignment

The audited defaults are BLOSUM62, gap opening `-10.0`, gap extension `-0.5`, and `max_gap=50`. The residue-alignment routine uses backward dynamic programming, chooses the best entry point over the score matrix, and permits alignment termination. Internal gap scoring is of the form

$$
\text{gap score}(L)=\text{gap}+\text{extend}\,(L-1).
$$

This is best described as a **local, Smith–Waterman-like residue alignment with PyMOL-specific gap/search controls**. It should not be treated as an unconditional promise of equivalence to an unrestricted textbook affine-gap implementation: `max_gap`, `max_skip`, and the structural-mode options affect the permitted paths and scores.[^wrapper-align][^match]

The critical practical point is that coordinate residuals are not fed back into a fresh sequence alignment during the later rejection loop. The residue alignment supplies the initial atom-pair pool once.[^executive]

### 2.2 Residue pairs become atom pairs

`ExecutiveAlign` passes the residue pairs to `SelectorCreateAlignments`. The latter finds corresponding selected atoms within those residue pairs; the call does not require identical residue names. Thus, for example, backbone atoms can be paired across a substitution, while atoms without a suitable selected counterpart are absent from the fitting list.[^executive][^selector]

A whole-protein selection therefore does **not** mean an implicit Cα-only fit. To make each observed residue contribute one coordinate, explicitly select Cα atoms. A backbone or all-atom selection has a different weighting interpretation because residues can contribute different numbers of paired atoms.[^selector][^wrapper-align]

For a reproducible comparison, preprocess alternate conformations, choose the intended chains and domains, and establish a missing-coordinate policy. Equal atom counts alone do not prove correspondence. These are workflow recommendations motivated by the way the atom-pair lists are constructed, not additional guarantees supplied by `align`.

### 2.3 What changes for `super`

`super` feeds the same overall C++ alignment machinery with different initial scoring parameters. Its defaults include `seq=0.0`, structural neighborhood parameters, and a local structural window; the wrapper also uses different gap penalties. The implementation can combine sequence and structural factors when these parameters are changed.[^wrapper-super][^executive][^match]

Consequently, `super` is not merely `align` with additional rejection cycles. It can start from a different residue correspondence. **Once its atom pairs reach `ExecutiveRMS`, however, the normalized-distance rejection mechanism is shared.**[^executive][^rejection]

## 3. The exact outer refinement rule

The decisive code is in `ExecutiveRMS`, in the audited file's approximately 10995–11065 region. For each current mobile coordinate, PyMOL applies the current trial transformation, calculates its distance from the corresponding target coordinate, and tests:

```cpp
if ((diff3f(v1, v2) / rms) > refine)
```

Here `refine` receives the command's `cutoff` parameter. The current RMSD is the RMSD calculated after fitting the current pair set.[^rejection]

Let $S_k$ be that set, $(R_k,t_k)$ its fitted transformation, and

$$
r_k=\sqrt{\frac{1}{|S_k|}\sum_{i\in S_k}d_{i,k}^{2}},
\qquad d_{i,k}=\|R_kx_i+t_k-y_i\|.
$$

The update is

$$
S_{k+1}=\{i\in S_k:d_{i,k}\le c\,r_k\},
$$

followed by a new fit on $S_{k+1}$. This equation is a direct translation of the inspected code, not a claim about an undocumented statistical model.[^rejection]

### 3.1 Pseudocode reflecting the implementation

```text
pairs = initial matched atom pairs
transform, rmsd = fit(pairs)
record initial pair count and initial fitted RMSD

for pass_number in 1 ... cycles:
    if cutoff or rmsd is numerically tiny:
        stop

    residuals = distances under the current transform
    survivors = current pairs with residual <= cutoff * rmsd
    # All decisions in this pass use the same transform and RMSD.

    if survivors is empty:
        fail: no atoms left after refinement

    had_rejections = survivors differs from current pairs
    pairs = survivors
    transform, rmsd = fit(pairs)
    record final statistics and this pass number

    if not had_rejections:
        stop

apply final transform if requested
return final statistics
```

The source also compacts the atom-identity arrays along with the coordinates, allowing an alignment object to record the surviving pairs.[^rejection][^output]

### 3.2 Consequences that are easy to miss

**Rejection is atom-wise.** If a residue contributes N, Cα, C, O, and side-chain atoms, some pairs can be rejected while others from the same residue remain. Selecting Cα makes the process effectively residue-wise only because there is normally one selected coordinate per observed residue.[^selector][^rejection]

**Rejected pairs are permanently removed during this call.** The routine compacts the coordinate arrays and operates on their reduced length. It neither reintroduces previously removed pairs nor substitutes nearby atoms.[^rejection]

**The test is relative, not absolute.** If the current RMSD is 4 Å, `cutoff=2` initially tolerates residuals up to 8 Å. After a later fit reaches 0.5 Å RMSD, the same parameter corresponds to a 1 Å threshold. These numbers are illustrative consequences of the source expression.

**It is not a conventional z-score.** The denominator is an RMS residual magnitude, not the sample standard deviation of distances about their mean. Calling the process “two-sigma clipping” without defining that distinction would be misleading.

**A final inspection pass can count as a cycle without removing atoms.** When refinement is active, the code can refit and update `n_cycles_run` even when all pairs survive that pass. Conversely, `cycles=5` is a maximum, not a promise of five rounds of removal.[^rejection]

**The cycle cap is not a fixed-point guarantee.** The final refit after the last allowed removal can change residuals and the RMSD threshold; another inspection could then remove additional pairs. There is no outer-loop test based on a minimum RMSD improvement.[^rejection]

**The printed rejection-line RMSD is the preceding fit's RMSD.** That message is emitted before the fit on the newly reduced pair set. The final reported statistics are updated afterward.[^rejection]

**The default fitting weights here are equal per atom.** `ExecutiveRMS` passes a null weight pointer to the fitting routine. This path does not automatically give each residue equal weight or weight atoms by mass, occupancy, B-factor, or prediction confidence.[^rejection][^matrix]

### 3.3 A numerical illustration, independently executed

I constructed 20 synthetic, symmetrically arranged coordinate pairs. Sixteen pairs had radial residuals of 0.1 units, two had residuals of 1 unit, and two had residuals of 6 units. Symmetry makes the identity rotation and zero translation the least-squares solution. Using the inspected relative threshold with `cutoff=2` gives:

| Inspection | Pairs entering inspection | Fitted RMSD | Threshold, `2 × RMSD` | Pairs removed |
|---|---:|---:|---:|---:|
| 1 | 20 | 1.925617 | 3.851234 | 2 |
| 2 | 18 | 0.346410 | 0.692820 | 2 |
| 3 | 16 | 0.100000 | 0.200000 | 0 |

The second pass removes the 1-unit residuals even though an absolute 2-unit cutoff would retain them. This illustrates the distinction without relying on any protein-specific assumptions. The independent NumPy calculation reproduced the displayed values; it does not prove that a particular installed PyMOL binary has the same implementation or numerical tolerances.

## 4. Interpreting the returned RMSD and counts

The Python result for `align` and `super` contains seven entries:[^binding]

| Index | Meaning | Important qualification |
|---:|---|---|
| 0 | Final RMSD | Fitted RMSD on surviving atom pairs |
| 1 | Final atom-pair count | Not necessarily a residue count |
| 2 | Refinement cycles run | May include a pass with no rejection |
| 3 | Initial RMSD | Already fitted, but before atom rejection |
| 4 | Initial atom-pair count | After correspondence construction |
| 5 | Raw residue-alignment score | Not an RMSD or a normalized similarity probability |
| 6 | Number of aligned residues | Recorded from the initial residue alignment, not recomputed as the final surviving-core residue count |

The initial RMSD therefore should not be described as the RMSD of the two original, untransformed coordinate frames. To measure the latter, use a known correspondence with a no-fit evaluation.[^rejection][^wrapper-other]

An alignment object created by the refined call contains the surviving atom pairs. Preserve a separate `cycles=0` alignment object when the initial correspondence is needed for later evaluation.[^output]

### Four useful quantities to report separately

My recommended reporting vocabulary is:

| Quantity | Fitting set | Evaluation set | Scientific question |
|---|---|---|---|
| Untrimmed fitted RMSD | All predefined paired atoms | Same full pair set | Overall least-squares agreement |
| Trimmed/core RMSD | Surviving core | Surviving core | How closely a selected rigid core agrees |
| Full-pair RMSD under core fit | Surviving core | Original full pair set | How much excluded regions deviate from the core's reference frame |
| Region RMSD under framework fit | Predefined framework | Predefined loop or interface | How a region changes relative to a common scaffold |

The last two **must not perform a second fit on the evaluation region**. Otherwise the question changes. This is especially relevant for antibody CDRs, flexible termini, interdomain orientations, and binder interfaces.

## 5. Other commands are not interchangeable

| Command | Correspondence / fitting behavior | Main distinction |
|---|---|---|
| `align` | Sequence-derived correspondence, fit, optional rejection | Defaults to five rejection cycles |
| `super` | Structurally informed initial residue correspondence, same fit/rejection backend | Different initial mapping from `align` |
| `fit` | Uses matching atoms or an explicit matching policy | Defaults to zero rejection cycles |
| `rms` | Calculates a trial fitted RMSD without moving the model | Not a current-pose RMSD |
| `rms_cur` | Calculates RMSD without fitting | Requires a trustworthy correspondence already |
| `pair_fit` | Fits explicitly supplied atom pairs / ordered selections | Useful when correspondence is externally defined |
| `cealign` | Separate combinatorial-extension structure alignment | Default guide atoms; argument order is **target, mobile** |
| `usalign` at the audited revision | Separate TM-score-oriented structural alignment | Availability must be checked in the installed version |

The distinctions through `cealign` are explicit in the Python wrappers.[^wrapper-other][^wrapper-super] CE has its own published algorithm and is not an invocation of the `ExecutiveRMS` rejection loop with different defaults.[^ce]

The audited source also contains a `usalign` wrapper using protein Cα / nucleic-acid guide atoms. Its existence in this source snapshot does not establish availability in an older packaged PyMOL installation. The US-align paper describes a TM-score objective and heuristic correspondence search, which are conceptually different from reporting a trimmed least-squares RMSD.[^wrapper-super][^usalign]

## 6. Two further reproducibility traps

### 6.1 Numerical fit iterations are not rejection cycles

`MatrixFitRMSTTTf` consults numerical-fitting settings such as `fit_iterations` and `fit_tolerance`. These concern finding a transformation for a fixed coordinate set. The outer `cycles` parameter controls atom rejection and repeated calls to the fitter. They are different loops.[^matrix][^rejection]

The source also offers a `fit_kabsch` switch, whose default is zero in `SettingInfo.h`. An independent SVD implementation can reproduce the mathematical rigid-fit objective without reproducing all numerical behavior of PyMOL's default solver. Near-degenerate coordinate sets and threshold-boundary cases deserve explicit regression tests.[^matrix][^settings]

### 6.2 Multi-state objects require an explicit policy

The Python defaults `mobile_state=0` and `target_state=0` are forwarded as negative internal state indices. In the inspected `ExecutiveRMS` coordinate-collection path, negative states use coordinates accumulated across states and averaged per atom. That is **not the same as independently fitting each state pair and averaging the RMSDs**.[^wrapper-align][^states]

For a single structure comparison, explicitly pass state 1 or another intended state. For ensembles, iterate the desired state pairs and define whether every model uses its own fit or a shared reference-frame fit. Record object transformation settings as well: PyMOL can represent movement through coordinates or object/state matrices.[^states][^output]

## 7. A practical PyMOL audit pattern

The following example assumes two clean, single-chain objects named `mobile` and `target`, one intended conformer per atom, and state 1. It is an example for execution inside PyMOL; it was not run against a PyMOL binary in this research session.

```python
from pymol import cmd
import numpy as np

mobile_selection = "mobile and polymer.protein and name CA"
target_selection = "target and polymer.protein and name CA"

# Save the initial mapping without moving either object.
initial = cmd.align(
    mobile_selection, target_selection,
    cycles=0, transform=0, object="all_pairs", reset=1,
    mobile_state=1, target_state=1,
)
original_pairs = cmd.get_raw_alignment("all_pairs")

# Fit and reject using the same sequence-alignment settings.
refined = cmd.align(
    mobile_selection, target_selection,
    cutoff=2.0, cycles=5, transform=1,
    object="core_pairs", reset=1,
    mobile_state=1, target_state=1,
)

# Evaluate every originally paired CA under the final core transform.
# Atom indices must remain unchanged between mapping and evaluation.
squared_distances = []
for column in original_pairs:
    pair = dict(column)  # raw alignment entries are (object_name, atom_index)
    if "mobile" not in pair or "target" not in pair:
        continue
    x = np.asarray(cmd.get_atom_coords(
        f"mobile and index {pair['mobile']}", state=1
    ), dtype=float)
    y = np.asarray(cmd.get_atom_coords(
        f"target and index {pair['target']}", state=1
    ), dtype=float)
    squared_distances.append(float(np.dot(x - y, x - y)))

if not squared_distances:
    raise RuntimeError("No complete original atom pairs were recovered")

print("Initial fitted RMSD / pairs:", initial[0], initial[1])
print("Core fitted RMSD / pairs:", refined[0], refined[1])
print("Retained pair fraction:", refined[1] / refined[4])
print("Original-pair RMSD under core fit:",
      float(np.sqrt(np.mean(squared_distances))))
```

For a production test, also save the exact selections, PyMOL version/build, scheme-derived residue mapping where applicable, states, gap parameters, fitting settings, surviving pair IDs, and transformation. Add a synthetic scaling test: under the relative rule, uniform coordinate scaling should preserve rejection membership apart from numerical-tolerance effects. This is a recommendation for validating the installed implementation, not a claim that the binary test was performed here.

## 8. Recommended use in protein and antibody design

For a near-identical design and reference, establish atom correspondence explicitly, calculate an untrimmed Cα or backbone RMSD, and retain residue-wise residuals. For divergent proteins, choose sequence-based or structure-based correspondence deliberately rather than assuming the command with the lowest RMSD found the biologically correct alignment.

For antibody comparisons, define the framework and CDR evaluation sets using an explicit numbering scheme and CDR convention. Fit a predefined framework, evaluate the CDRs without refitting, and separately quantify VH–VL orientation changes when relevant. An automatically clipped whole-variable-domain RMSD may remove precisely the loop changes the design experiment is intended to measure.

For a rigid-core visualization, `align` or `super` with rejection is useful. For a scientific ranking, treat the resulting core selection as part of the result, not as an invisible preprocessing step.

## 9. Evidence limitations

I did not identify a dedicated peer-reviewed publication that specifies this exact PyMOL normalized-distance rejection implementation. The classical least-squares, sequence-alignment, CE, and US-align publications explain relevant algorithm families, but they should not be cited as proof of PyMOL's particular clipping threshold. The pinned C++ code is the decisive evidence for that behavior.

The audit covers the stated open-source revision. It does not establish that every historical release, commercial build, third-party patch, or future release behaves identically. A source–binary mismatch should be resolved by running a small regression fixture on the deployed executable.

## References and audit links

[^revision]: Schrödinger, **PyMOL open-source repository**, cutoff-constrained commit query: <https://api.github.com/repos/schrodinger/pymol-open-source/commits?until=2026-09-08T23:59:59Z&per_page=1>. Returned revision: <https://github.com/schrodinger/pymol-open-source/commit/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69>.
[^wrapper-align]: PyMOL, pinned `modules/pymol/fitting.py`, `align` defaults, docstring, and argument forwarding: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/modules/pymol/fitting.py#L370-L463>. In particular, the docstring calls the cutoff an Å value, whereas the rejection source below normalizes by RMSD.
[^rejection]: PyMOL, pinned `layer3/Executive.cpp`, initial fit and normalized-distance rejection: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/layer3/Executive.cpp#L10990-L11090>.
[^executive]: PyMOL, pinned `layer3/Executive.cpp`, especially `ExecutiveAlign` and its calls to `MatchAlign`, `SelectorCreateAlignments`, and `ExecutiveRMS`: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/layer3/Executive.cpp>.
[^selector]: PyMOL, pinned `layer3/Selector.cpp`, `SelectorCreateAlignments`: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/layer3/Selector.cpp>.
[^kabsch]: Kabsch W. **A solution for the best rotation to relate two sets of vectors.** *Acta Crystallographica A* 32, 922–923 (1976). DOI: <https://doi.org/10.1107/S0567739476001873>.
[^matrix]: PyMOL, pinned `layer0/Matrix.cpp`, `MatrixFitRMSTTTf`: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/layer0/Matrix.cpp>.
[^settings]: PyMOL, pinned `layer1/SettingInfo.h`, `fit_kabsch` default, and setting documentation: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/layer1/SettingInfo.h>; <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/data/setting_help.csv>.
[^match]: PyMOL, pinned `layer0/Match.cpp`, `MatchAlign`: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/layer0/Match.cpp#L370-L660>.
[^binding]: PyMOL, pinned `layer4/Cmd.cpp`, construction of the result using `final_rms`, `final_n_atom`, `n_cycles_run`, `initial_rms`, `initial_n_atom`, `raw_alignment_score`, and `n_residues_aligned`: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/layer4/Cmd.cpp>.
[^wrapper-super]: PyMOL, pinned `modules/pymol/fitting.py`, `cealign`, `usalign`, and `super`: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/modules/pymol/fitting.py#L27-L373>.
[^output]: PyMOL, pinned `layer3/Executive.cpp`, survivor alignment object and final transformation handling: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/layer3/Executive.cpp#L11075-L11230>.
[^wrapper-other]: PyMOL, pinned `modules/pymol/fitting.py`, `fit`, `rms`, `rms_cur`, and `pair_fit`: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/modules/pymol/fitting.py#L610-L825>.
[^ce]: Shindyalov IN, Bourne PE. **Protein structure alignment by incremental combinatorial extension (CE) of the optimal path.** *Protein Engineering* 11, 739–747 (1998). DOI: <https://doi.org/10.1093/protein/11.9.739>. PubMed: <https://pubmed.ncbi.nlm.nih.gov/9796821/>.
[^usalign]: Zhang C, Shine M, Pyle AM, Zhang Y. **US-align: universal structure alignments of proteins, nucleic acids, and macromolecular complexes.** *Nature Methods* 19, 1109–1115 (2022). DOI: <https://doi.org/10.1038/s41592-022-01585-1>.
[^states]: PyMOL, pinned `layer3/Executive.cpp`, `ExecutiveRMS` coordinate acquisition and state averaging: <https://github.com/schrodinger/pymol-open-source/blob/5e8bfca5a7f5dc4d5e7f84fa1d15af707cc86e69/layer3/Executive.cpp#L10580-L10730>.
