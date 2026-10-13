# Project notes: challenges, mistakes and lessons

Honest record of what went wrong and what was learned while building the pipeline.

## What was built

SigXTalk applied to GSE174367 (Morabito et al. 2021) with microglia as the receiver cell type, run separately for
Alzheimer's disease (AD) and Control nuclei, with the same target genes and thresholds in both groups.

## Challenges and how they were handled

| Problem | What happened | Fix / decision |
|---|---|---|
| Memory | CellChat on tens of thousands of nuclei exceeds R's default vector-memory cap on macOS | Cap lifted for the run (`RAISE_R_MEMORY_CAP`); every step saves a checkpoint to `work/` so an interrupted run resumes |
| Downsampling | Early runs capped each cell type at 1,500 nuclei | Final run uses all 61,472 nuclei; the cap changed which cells informed the model |
| Reproducibility | The random-forest step in the package function was not seeded | Wrote `PRS_calc_seeded()` that seeds each parallel worker (`PRS_SEED = 42`) |
| Few receptors survive | Control kept only one receptor (PTPRM) | Investigated: not a threshold problem (see below) |
| Missing object | A diagnostic failed with "object 'SeuratObj' not found" | The pipeline removes the full object after building per-group data to save memory; diagnostics need Sections 2-4 re-run |

## Why so few receptors survive (the main finding about the method)

SigXTalk keeps a ligand-receptor pair only if **both** genes are in the receiver's expression matrix. That matrix
contains only the ~2,000 variable genes that were scaled, further reduced to genes detected in at least 5% of
microglia (444 genes for Control). Genes such as APP, SORL1, SPP1, ITGB2 and NOTCH2 are not in it, so pairs that depend
on them cannot appear, and the PTPRC -> MRC1 weight in Control is already zero at the CellChat stage. Lowering the
ligand-receptor threshold therefore cannot add MRC1 to Control. The effective bottleneck is the gene universe. Scaling
more genes (all genes, or a larger variable set) is the most important next step. The SigXTalk paper does not report a
sensitivity analysis of its 5% / 10% thresholds, so no threshold in this project is claimed to be optimal.

## Mistakes made along the way and corrections

* The direction of the alluvial plot was first described the wrong way round; it reads receptor -> intermediate ->
  target, top to bottom.
* MRC1 was first discussed as a target; it is a receptor.
* The lost PTPRC -> MRC1 link in Control was first attributed to a threshold crossing; it is zero in CellChat.
* A claim that AD "loses" APP -> SORL1/CD74 signaling became provisional after the checkpoint mismatch below.
* Several statements about run time and about seed changes in `main.py` were removed from the documentation because
  they could not be verified.

## Unresolved issue

The user's Control CellChat checkpoint, read through `Extract_LR_Prob()`, has no APP rows, but the saved comparison
table lists APP -> SORL1 and APP -> CD74 for Control. PTPRM -> PTPRM matches. Suggested check:

```r
ad_lr <- Extract_LR_Prob(readRDS("work/<AD cellchat file>.rds"), "MG")
ad_lr[ad_lr$From == "APP", ]
```

If APP is present for AD but absent for Control in the checkpoints, the table was produced from a different Control
object than the one stored; re-run the comparison from the current checkpoints.

## Interpretation limits

* Pathways come from expression co-variation; they are hypotheses, not demonstrated signaling.
* Weights are unsigned.
* More routes in AD are partly expected because AD has twice as many microglia (2,742 vs 1,384) and the number of routes is
  capped by the number of surviving receptors.
* In the direction check, the genes that are higher in AD are detected in only a minority of microglia (about 8-35%),
  which suggests a shift in microglial state mix rather than a uniform up-regulation.
* DE p-values use nuclei as replicates; with donors as replicates they would be far less significant.

## What I would do differently

1. Scale all genes (or at least several thousand more) before SigXTalk, then compare the receptor list.
2. Run a threshold sensitivity analysis (expression filter, ligand-receptor cut-off, HGNN activity cut-off).
3. Use donor-level (pseudobulk) tests for the target genes.
4. Re-derive the ligand-receptor comparison table from the stored checkpoints in the same run as the figures.
5. Repeat for a second cell type to check that conclusions are not microglia-specific.

## Provenance note

The stored results were produced by the equivalent plain R script that this R Markdown file was converted from.
The R Markdown file contains the same code with repository-relative paths. It has been syntax-checked but should
be knitted once end-to-end to confirm it reproduces the stored files.
