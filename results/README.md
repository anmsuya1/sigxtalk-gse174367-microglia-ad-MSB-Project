# Results: guide to every file

Everything in this folder was produced by `analysis/01_sigxtalk_pipeline.Rmd` (and, for one table,
`analysis/02_de_direction_check.Rmd`). These are the files behind the numbers in the top-level `README.md`.
Re-knitting the pipeline overwrites them.

Vocabulary: **receptor** = the microglial receptor that receives a signal; **SSC** = signaling-and-sensing component
(the intermediate/transcription factor between receptor and target); **target** = AD-vs-Control differentially
expressed microglial gene; **weight / PRS** = pathway regulatory strength from the random forests (higher = the
route explains more of the target's variation; it is unsigned).

## Figures (`figures/GSE174367_<group>_full/`)

One folder per group (`Control`, `AD`). The file names are identical so the two groups can be compared side by side.

| File | What it shows | How to read it |
|---|---|---|
| `12a_CCI_overview_circle.png` | Cell-cell communication between all seven cell types (CellChat) | Arrow thickness = summed communication strength |
| `12a_LR_circle_MG.png` | Ligand-receptor pairs that send signals **into microglia** | Which sender cell type uses which ligand-receptor pair to reach MG |
| `12b_crosstalk_counts_distribution.png` | How many routes reach each target gene | Distribution over targets |
| `12b_crosstalk_counts_topgenes.png` | The targets reached by the most routes | A target with many routes is regulated by many receptor/intermediate combinations |
| `12c_receptor_target_heatmap.png` | Receptor-by-target strength (**AD only**) | Which receptor drives which target |
| `12d_alluvial_top3.png` | Strongest routes as flows: receptor -> intermediate -> target | Flow width = route weight; read top to bottom |
| `12d_heatmap_DOCK8.png` | Receptor-by-intermediate strengths for the target DOCK8 (**AD only**) | Which intermediates carry each receptor's signal to DOCK8 |
| `12d_fidelity_specificity_DOCK8.png` (AD), `12d_fidelity_specificity_BLNK.png` (Control) | Fidelity (share of the target's regulation passing through a route) and specificity (how well the signal avoids non-targets) | Routes high on both are the most selective |

**Why Control has no `12c` and no `12d_heatmap`:** both need at least two receptors, and only one receptor (PTPRM)
survives in Control. This is a result, not an error. Both groups use the same targets and thresholds.

## Tables (`tables/`)

| File | Content |
|---|---|
| `results_Control_full.csv`, `results_AD_full.csv` | Final pathways after all filters: `Receptor`, `SSC`, `Target`, `Weight` |
| `pathways_merged_AD_vs_Control.csv` | Both groups side by side (`Weight_AD`, `Weight_Control`; `NA` = route absent in that group) |
| `pathways_AD_specific.csv`, `pathways_Control_specific.csv` | Routes present in only one group |
| `pathway_counts_per_target_AD.csv`, `pathway_counts_per_target_Control.csv` | Number of surviving routes per target gene |
| `CCI_celltype_weights_AD_vs_Control.csv` | Summed CellChat strength for each sender-receiver cell-type pair, and AD minus Control |
| `LR_pairs_into_MG_AD_vs_Control.csv` | Ligand-receptor weights into microglia in each group. **The Control column is provisional**, see below |
| `pathway_genes_DE_direction_MG.csv` | Direction of change (AD vs Control) of every gene in the final pathways, in microglia only. Wilcoxon test on single nuclei, so the p-values are optimistic (nuclei are not independent donors) |

## Raw HGNN output (`hgnn_outputs/`)

`pathways_MG.csv` in each group folder lists every candidate Receptor -> TF -> target link that was given to the
neural network, with its predicted activity (`pred_label`). Only links with `pred_label > 0.75` continue to the
random-forest step (11 in Control of 699 candidates, 58 in AD of 1,404).

## `sessionInfo_full_run.txt`

R version and package versions of the run that produced these files.

## Known issue (provisional table)

When spot-checking, the Control CellChat checkpoint did not reproduce the APP -> SORL1 and APP -> CD74 rows that
appear in `LR_pairs_into_MG_AD_vs_Control.csv`, whereas PTPRM -> PTPRM did match. Until the Control side of this table
is re-derived, do not draw AD-versus-Control conclusions from the APP-based rows. The pathway tables, figures and
cell-type communication table do not depend on this table.
