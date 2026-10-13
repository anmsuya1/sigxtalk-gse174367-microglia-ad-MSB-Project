## ============================================================================
## SigXTalk on GSE174367 (Morabito et al. 2021, Nat Genet)
## Single-nucleus RNA-seq, human prefrontal cortex, AD vs Control
##
## FULL-CELL VERSION: NO DOWNSAMPLING.
## Every cell of every type goes into CellChat / SigXTalk. There is no cell cap,
## no Kennard-Stone step, and no random sampling anywhere in this script.



## --- Settings you can flip --------------------------------------------------
target_type        <- "MG"    # receiver cell type: MG, ASC, ODC, EX, INH, OPC, PER.END
USE_CHECKPOINTS    <- TRUE    # reuse finished conditions found on disk. To force a
                              # full recompute, delete the *.rds files in work_dir.
RAISE_R_MEMORY_CAP <- TRUE    # macOS R caps vector memory at 16 GB by default.
                              # TRUE lifts the cap so macOS can swap to disk instead
                              # of R stopping with "vector memory limit reached".
                              # Slower once it swaps, but it can finish. FALSE = keep cap.
DROP_COUNTS_LAYER  <- FALSE    # raw counts are not used after normalization; dropping
                              # that layer frees memory. Skipped automatically if it fails.
PRS_SEED           <- 42

## --- 0. Set up a clean project folder automatically ----------------------
downloads_dir <- path.expand("~/Downloads")
work_dir      <- path.expand("~/SigXTalk_GSE174367_FULL")
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

h5_file   <- "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"
meta_file <- "GSE174367_snRNA-seq_cell_meta.csv"

## Copy the two GEO files from Downloads into the project folder (only if
## not already copied)
if (!file.exists(file.path(work_dir, h5_file))) {
  ok <- file.copy(file.path(downloads_dir, h5_file), work_dir)
  if (!ok) stop("Couldn't find ", h5_file, " in ", downloads_dir,
                ". Make sure it's downloaded there first.")
}
if (!file.exists(file.path(work_dir, meta_file))) {
  ok <- file.copy(file.path(downloads_dir, meta_file), work_dir)
  if (!ok) stop("Couldn't find ", meta_file, " in ", downloads_dir,
                ". Make sure it's downloaded there first.")
}

## Fetch main.py straight from the SigXTalk repo 
if (!file.exists(file.path(work_dir, "main.py"))) {
  download.file(
    "https://raw.githubusercontent.com/LithiumHou/SigXTalk/master/vignettes/main.py",
    destfile = file.path(work_dir, "main.py")
  )
}

stopifnot(file.exists(file.path(work_dir, h5_file)),
          file.exists(file.path(work_dir, meta_file)),
          file.exists(file.path(work_dir, "main.py")))
cat("Working directory ready:", work_dir, "\n")
print(list.files(work_dir))

## Memory report (and optionally lift R's vector-memory cap)
phys_gb <- tryCatch(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 1024^3,
                    error = function(e) NA_real_, warning = function(w) NA_real_)
cat("\nPhysical RAM on this machine: ",
    ifelse(is.na(phys_gb), "unknown", paste0(round(phys_gb, 1), " GB")), "\n", sep = "")
cat("R vector-memory cap before change (MB):", tryCatch(mem.maxVSize(), error = function(e) NA), "\n")
if (RAISE_R_MEMORY_CAP) {
  try(mem.maxVSize(Inf), silent = TRUE)
  cat("R vector-memory cap after change  (MB):", tryCatch(mem.maxVSize(), error = function(e) NA), "\n")
}

## --- 1. Load libraries -----------------------------------------------------
## (install first if needed: install.packages("Seurat"); devtools::install_github("sqjin/CellChat");
##  devtools::install_github("LithiumHou/SigXTalk"))
library(Seurat)
library(dplyr)
library(CellChat)
library(SigXTalkR)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(foreach)   

## --- 2. Load the raw data and PEEK
mat <- Read10X_h5(h5_file)
cat("\nRaw matrix dimensions (genes x cells):", nrow(mat), "x", ncol(mat), "\n")
cat("First 5 genes x first 5 cells (raw UMI counts):\n")
print(mat[1:5, 1:5])

meta <- read.csv(meta_file)
cat("\nMetadata dimensions (cells x columns):", nrow(meta), "x", ncol(meta), "\n")
cat("First 5 rows of metadata:\n")
print(head(meta, 5))

## --- 3. Sanity check + fix the barcode mismatch ----------------------------
## The matrix has 61,770 barcodes; the metadata only annotates 61,472 of them
## (298 cells in the matrix were excluded from the paper's final cell typing).

cat("\n--- Sanity check: matrix vs metadata ---\n")
cat("Barcodes in matrix:  ", ncol(mat), "\n")
cat("Barcodes in metadata:", nrow(meta), "\n")
cat("Matrix barcodes with NO metadata annotation:", sum(!colnames(mat) %in% meta$Barcode), "\n")
cat("Metadata barcodes missing from matrix (should be 0):", sum(!meta$Barcode %in% colnames(mat)), "\n")

mat <- mat[, colnames(mat) %in% meta$Barcode]     # drop the 298 unannotated cells
rownames(meta) <- meta$Barcode
meta <- meta[colnames(mat), ]                      # align order to the matrix

stopifnot(all(colnames(mat) == rownames(meta)))
cat("\nAfter alignment - matrix:", ncol(mat), "cells | metadata:", nrow(meta), "cells (should match)\n")

## --- 4. Build the Seurat object -------------------------------------------
## No clustering step needed - Cell.Type is already provided by the authors.
SeuratObj <- CreateSeuratObject(counts = mat, meta.data = meta,
                                 min.cells = 3, min.features = 200)
rm(mat); gc()    # the raw matrix is now duplicated inside SeuratObj - free it right away
cat("\nSeurat object after QC filtering:", ncol(SeuratObj), "cells x", nrow(SeuratObj), "genes\n")


SeuratObj <- NormalizeData(SeuratObj)
SeuratObj <- FindVariableFeatures(SeuratObj)
SeuratObj <- ScaleData(SeuratObj, features = VariableFeatures(SeuratObj))
gc()

Idents(SeuratObj) <- SeuratObj$Cell.Type       # ODC / EX / INH / ASC / MG / OPC / PER.END
cat("\nCell type x Diagnosis breakdown (ALL of these cells are used - no downsampling):\n")
print(table(SeuratObj$Cell.Type, SeuratObj$Diagnosis))

## Diagnostic: Prepare_Input reads scale.data, so the genes in scale.data are the
## only genes SigXTalk can use as receptors / TFs / targets. This prints how
## many there are and whether the genes from earlier results are present.
sd_genes <- tryCatch(rownames(GetAssayData(SeuratObj, assay = "RNA", layer = "scale.data")),
                     error = function(e) NULL)
if (!is.null(sd_genes)) {
  cat("\nGenes in scale.data (the gene universe SigXTalk's Prepare_Input can see):", length(sd_genes), "\n")
  key_genes <- c("CSF1R", "TREM2", "SORL1", "CD74", "APP", "MRC1", "SYK", "ESR1", "NRP1", "NOTCH2", "CNTN1", "PALD1")
  print(data.frame(gene = key_genes, in_scale_data = key_genes %in% sd_genes))
}

## --- 5. Pick the receiver cell type ----------------------------------------
## (target_type is set in the Settings block at the top.) Condition-specific subsets are built in
## Section 9, right before each condition runs, to keep peak memory down.
cat("\nReceiver cell type:", target_type, "\n")

## --- 6. Define target genes = AD-vs-Control DE genes WITHIN the receiver

mg_only <- subset(SeuratObj, subset = Cell.Type == target_type)
Idents(mg_only) <- mg_only$Diagnosis
AD_DEGs <- FindMarkers(mg_only, ident.1 = "AD", ident.2 = "Control",
                        min.pct = 0.1, logfc.threshold = 0.25)
TG_used <- AD_DEGs %>% dplyr::filter(p_val_adj < 0.05) %>% rownames()
cat("\nNumber of AD-vs-Control DE genes found in", target_type, ":", length(TG_used), "\n")
rm(mg_only, AD_DEGs); gc()

## --- 7. Load SigXTalk's built-in human prior networks ----------------------
data("RTF_human")   # Receptor -> TF
data("TFT_human")   # TF -> Target gene

## Only gene NAMES are needed here. rownames(SeuratObj) reads them from the
## object's feature list without touching the expression matrix earlier versions 
## materialized the whole matrix and hit the 16 GB vector-memory limit).
all_genes <- rownames(SeuratObj)
RecTFDB <- RTF_human %>% distinct(from, to, .keep_all = TRUE) %>%
  Filter_DB(all_genes)
TFTGDB  <- TFT_human %>% distinct(from, to, .keep_all = TRUE) %>%
  Filter_DB(all_genes)

## --- 7b. Reproducibility patch for PRS_calc() -------------------------------
## PRS_calc() fits one random forest (ranger) per TF-TG pair on its OWN internal
## parallel cluster (snow::makeSOCKcluster), and nothing seeds those workers, so
## identical inputs can give different PRS weights run to run. This is an
## identical copy of PRS_calc()'s source with one added line,
## parallel::clusterSetRNGStream(cl, seed), right after the cluster is created.
## The unexported SigXTalkR helpers (Run_Regression,
## Calculate_Importance) are called with ::: .
PRS_calc_seeded <- function(Exp_clu, RTT_results, cutoff = 0.15, seed = 42) {

  message("Preprocessing data...\n")
  colnames(RTT_results)[1:3] <- c("Receptor", "TF", "TG")
  RecAll <- RTT_results$Receptor %>% unique()
  TFAll <- RTT_results$TF %>% unique()
  TGAll <- RTT_results$TG %>% unique()

  RTT_results$TFTG <- paste0(RTT_results$TF, "_", RTT_results$TG)

  RecofTFTG <- split(RTT_results, RTT_results$TFTG)
  TFTG_All <- names(RecofTFTG)

  Ncores <- parallel::detectCores()
  cl <- snow::makeSOCKcluster(min(24, Ncores))
  parallel::clusterSetRNGStream(cl, seed)   ## <-- the actual fix
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(min = 0, max = length(TFTG_All), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  message(paste0("\n Estimating the nonlinear model, including ", length(TFTG_All), " TF-TG pairs. \n"))

  ModelList <- foreach::foreach(
    i = 1:length(TFTG_All), .packages = c(
      "dplyr", "ranger",
      "tidymodels", "rlang", "parsnip",
      "recipes"
    ),
    .export = c(
      "TFTG_All", "Exp_clu", "cutoff"
    ),
    .options.snow = opts, .errorhandling = "pass"
  ) %dopar% {

    tftg_name <- TFTG_All[i]
    tf_name <- strsplit(tftg_name, "_", fixed = TRUE)[[1]][1]
    tg_name <- strsplit(tftg_name, "_", fixed = TRUE)[[1]][2]
    message("Current TF-TG pair: ", tftg_name, "\n")

    Exp_tftg <- Exp_clu[tf_name, ] * Exp_clu[tg_name, ]
    Exp_tftg <- as.data.frame(Exp_tftg)
    tar_allrec <- RecofTFTG[[tftg_name]]$Receptor %>% unique()
    tar_allrec <- setdiff(tar_allrec, c(tf_name, tg_name))
    ncells <- sum(Exp_clu[tg_name, ] > 0)

    if (rlang::is_empty(tar_allrec)) {
      cat("No corresponding Rec, skip this TF... \n")
      next
    }
    if (sum(Exp_tftg > 1e-05) < cutoff * ncells) {
      cat("No enough expression, skip this TF... \n")
      next
    }

    Exp_rec_used <- Exp_clu[tar_allrec, ] %>% as.matrix()

    if (min(dim(Exp_rec_used)) > 1) {
      Exp_rec_used <- t(Exp_rec_used)
      Exp_rec_used <- scale(Exp_rec_used)
    }
    if (!max(Exp_rec_used)) {
      cat("No expression of receptors, skip this TF... \n")
      next
    }

    data_use <- cbind(Exp_tftg, Exp_rec_used)
    colnames(data_use)[1:2] <- c("y", tar_allrec[1])
    colnames(data_use) <- gsub("-", "_", colnames(data_use))
    SigXTalkR:::Run_Regression(data_use = data_use)
  }
  names(ModelList) <- TFTG_All

  close(pb)
  parallel::stopCluster(cl)
  rm(cl)
  gc()

  message("Estimating importance scores of Recs...\n")

  coefList <- SigXTalkR:::Calculate_Importance(ModelList)
  RecTFTG <- do.call(rbind, coefList) %>% as.data.frame()
  RecTFTG_final <- RecTFTG %>% tidyr::separate(To, into = c("SSC", "Target"), sep = "_")

  RecTFTG_final <- RecTFTG_final[, c("From", "SSC", "Target", "Weight")]
  colnames(RecTFTG_final) <- c("Receptor", "SSC", "Target", "Weight")

  return(RecTFTG_final)
}

## --- 8. Helper: run the full SigXTalk pipeline on one condition -----------
## Checkpoints EVERY intermediate object it computes 
run_sigxtalk <- function(obj, target_type, TGs, project_name, work_dir,
                          conda_env = "SigXTalk_py") {

  cell_anno <- data.frame(cell = colnames(obj), cluster = as.character(Idents(obj)))

  ## 8a. Infer cell-cell communication with CellChat (RAM-heavy step; it
  ##     sees every cell of every type here)
  lr_file <- file.path(work_dir, paste0("LR_original_", project_name, ".rds"))
  if (USE_CHECKPOINTS && file.exists(lr_file)) {
    cat("\n[", project_name, "] reusing saved CellChat result: ", lr_file, "\n", sep = "")
    LR_original <- readRDS(lr_file)
  } else {
    LR_original <- Infer_CCI(obj, cell_anno, use_spatial = FALSE, db_use = "human")
    saveRDS(LR_original, lr_file)
  }
  gc()

  ## 8b. Write the HGNN module's input files (inputs/ExpressionCount.csv,
  ##     inputs/RecTFDB.csv, inputs/TFTGDB.csv)
  input_dir <- file.path(work_dir, "inputs")
  dir.create(input_dir, showWarnings = FALSE)
  Prepare_Input(obj, target_type, TGs = TGs, CCC_results = LR_original,
                RecTFDB, TFTGDB, data_dir = input_dir,
                assay = "RNA", datatype = "scale.data",
                exp_threshold = 0.05, CCC_threshold = 0.05)

  ## 8c. Train the hypergraph neural network (Python side) to score which
  ##     Receptor -> TF -> Target paths are actually "activated"
  args <- c("--project", project_name, "--target_type", target_type,
            "--device", "cpu")           # forced to CPU on Apple Silicon regardless (no CUDA)
  Run_py_script("main.py", conda_env, args)

  ## 8d. Compute the Pathway Regulatory Strength (PRS) and filter
  output_dir <- file.path(work_dir, "outputs", project_name)
  RTFTG_results <- read.csv(file.path(output_dir, paste0("pathways_", target_type, ".csv")))
  RTFTG_results <- RTFTG_results[RTFTG_results$pred_label > 0.75, 1:3]
  Exp_clu <- Get_Exp_Clu(obj, clusterID = target_type, assay = "RNA",
                          datatype = "data", cutoff = 0.1)
  saveRDS(Exp_clu, file.path(work_dir, paste0("Exp_clu_", project_name, ".rds")))
  cat("\n[", project_name, "] ", target_type, " cells used for PRS: ", ncol(Exp_clu),
      " | genes: ", nrow(Exp_clu), " | HGNN-activated pathways (pred_label > 0.75): ",
      nrow(RTFTG_results), "\n", sep = "")

  ress <- PRS_calc_seeded(Exp_clu, RTFTG_results, cutoff = 0.1, seed = PRS_SEED)

  ## Checkpoint the raw PRS results BEFORE filtering, so a bug or argument
  ## mismatch in Filter_results() never costs re-running CellChat + HGNN
  ## training again - just reload this file and retry the filter step.
  saveRDS(ress, file.path(work_dir, paste0("ress_", project_name, ".rds")))

  results <- Filter_results(ress, PRS_thres = 0.01)   
  saveRDS(results, file.path(work_dir, paste0("results_", project_name, ".rds")))
  results
}

## --- 9. Run it once per condition, same target genes held constant --------
## Peak memory is kept down by: (1) dropping the unused raw-counts layer,
## (2) building the two condition subsets, then deleting the full object, and
## (3) deleting each condition's object as soon as its run has finished.
proj       <- c(Control = "GSE174367_Control_full", AD = "GSE174367_AD_full")
diag_value <- c(Control = "Control",                AD = "AD")

ckpt_files <- function(project_name) {
  file.path(work_dir, c(paste0("LR_original_", project_name, ".rds"),
                        paste0("Exp_clu_",     project_name, ".rds"),
                        paste0("ress_",        project_name, ".rds"),
                        paste0("results_",     project_name, ".rds")))
}
ckpt_complete <- function(project_name) USE_CHECKPOINTS && all(file.exists(ckpt_files(project_name)))

drop_counts <- function(o) {
  tryCatch({
    LayerData(o, assay = "RNA", layer = "counts") <- NULL
    cat("Dropped the raw counts layer (not used downstream) to save memory.\n")
    o
  }, error = function(e) {
    cat("Could not drop the counts layer (", conditionMessage(e),
        ") - continuing without; this only costs memory.\n", sep = "")
    o
  })
}

make_condition_obj <- function(full_obj, diagnosis_value) {
  o <- subset(full_obj, subset = Diagnosis == diagnosis_value)
  Idents(o) <- o$Cell.Type
  o
}

if (DROP_COUNTS_LAYER) SeuratObj <- drop_counts(SeuratObj)

cond_objs <- list()
for (k in names(proj)) {
  if (ckpt_complete(proj[[k]])) {
    cat("\n[", k, "] complete checkpoint found on disk - will load it instead of re-running.\n", sep = "")
  } else {
    cond_objs[[k]] <- make_condition_obj(SeuratObj, diag_value[[k]])
  }
}
rm(SeuratObj); gc()    # the full object is no longer needed - everything below uses the subsets

results <- list()
for (k in c("Control", "AD")) {        # Control first - it is the smaller condition
  res_file <- file.path(work_dir, paste0("results_", proj[[k]], ".rds"))
  if (ckpt_complete(proj[[k]])) {
    results[[k]] <- readRDS(res_file)
    cat("\n[", k, "] loaded from checkpoint: ", nrow(results[[k]]), " pathways\n", sep = "")
  } else {
    cat("\n================ Running ", k, " (all cells, no downsampling) ================\n", sep = "")
    cat("Cells per type in ", k, ":\n", sep = "")
    print(table(cond_objs[[k]]$Cell.Type))
    results[[k]] <- tryCatch(
      run_sigxtalk(cond_objs[[k]], target_type, TG_used, proj[[k]], work_dir),
      error = function(e) {
        cat("\n*** The ", k, " run stopped with: ", conditionMessage(e), "\n",
            "If this is a 'vector memory limit' / out-of-memory error, the full-cell\n",
            "run does not fit on this machine. Anything already finished is saved\n",
            "in ", work_dir, " (see the *.rds files).\n", sep = "")
        stop(e)
      }
    )
    cond_objs[[k]] <- NULL
    gc()
  }
}
results_Ctrl <- results[["Control"]]
results_AD   <- results[["AD"]]

## --- 10. Compare crosstalk rewiring between AD and Control -----------------
## Pathways with high PRS in AD but absent/low in Control (or vice versa) are
## candidate disease-specific signaling rewiring events.
##
## PRS_calc()'s actual output columns are Receptor, SSC, Target, Weight 
##"SSC" is the regulator column and "Weight" is the PRS score.
table_dir <- file.path(work_dir, "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n--- Pathway counts (all cells, no downsampling) ---\n")
cat("AD:      ", nrow(results_AD),   "pathways,", dplyr::n_distinct(results_AD$Target),   "distinct target genes\n")
cat("Control: ", nrow(results_Ctrl), "pathways,", dplyr::n_distinct(results_Ctrl$Target), "distinct target genes\n")

## per-gene diversity - how many genes have >=2 Receptors AND >=2 SSCs?
diversity_check <- function(results) {
  results %>%
    dplyr::group_by(Target) %>%
    dplyr::summarise(n_receptors = dplyr::n_distinct(Receptor),
                     n_ssc       = dplyr::n_distinct(SSC),
                     n_pathways  = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_pathways))
}
cat("\nPer-gene diversity - AD:\n");      print(diversity_check(results_AD))
cat("\nPer-gene diversity - Control:\n"); print(diversity_check(results_Ctrl))

ad_paths   <- results_AD   %>% dplyr::select(Receptor, SSC, Target, Weight)
ctrl_paths <- results_Ctrl %>% dplyr::select(Receptor, SSC, Target, Weight)

merged <- full_join(ad_paths, ctrl_paths,
                     by = c("Receptor", "SSC", "Target"),
                     suffix = c("_AD", "_Control"))
ad_specific   <- merged %>% dplyr::filter(!is.na(Weight_AD) & is.na(Weight_Control))
ctrl_specific <- merged %>% dplyr::filter(is.na(Weight_AD) & !is.na(Weight_Control))
shared        <- merged %>% dplyr::filter(!is.na(Weight_AD) & !is.na(Weight_Control))

cat("\nAD-specific crosstalk pathways:     ", nrow(ad_specific), "\n")
cat("Control-specific crosstalk pathways:", nrow(ctrl_specific), "\n")
cat("Shared (found in both):             ", nrow(shared), "\n")

## Which TARGET GENES are shared vs condition-specific
tg_ad <- unique(results_AD$Target); tg_ct <- unique(results_Ctrl$Target)
cat("\nTarget genes in both conditions:  ", paste(intersect(tg_ad, tg_ct), collapse = ", "), "\n")
cat("Target genes only in AD:          ", paste(setdiff(tg_ad, tg_ct), collapse = ", "), "\n")
cat("Target genes only in Control:     ", paste(setdiff(tg_ct, tg_ad), collapse = ", "), "\n")

## Consistency check between the pathway table and the alluvial plots: every
## pathway feeding the top genes, side by side for the two conditions.
top_genes <- unique(c(head(diversity_check(results_AD)$Target, 3),
                      head(diversity_check(results_Ctrl)$Target, 3)))
cat("\nAll pathways (both conditions side by side) for the top genes:\n")
print(merged %>% dplyr::filter(Target %in% top_genes) %>% dplyr::arrange(Target, Receptor, SSC))

write.csv(results_AD,     file.path(table_dir, "results_AD_full.csv"),      row.names = FALSE)
write.csv(results_Ctrl,   file.path(table_dir, "results_Control_full.csv"), row.names = FALSE)
write.csv(merged,         file.path(table_dir, "pathways_merged_AD_vs_Control.csv"), row.names = FALSE)
write.csv(ad_specific,    file.path(table_dir, "pathways_AD_specific.csv"),      row.names = FALSE)
write.csv(ctrl_specific,  file.path(table_dir, "pathways_Control_specific.csv"), row.names = FALSE)

## --- 10b. Numeric comparison of the CellChat results -----
## CellChat's edge widths are scaled within each plot, so two plots can look
## alike while the numbers differ. These tables compare the raw weights.
load_condition <- function(project_name) {
  list(LR  = readRDS(file.path(work_dir, paste0("LR_original_", project_name, ".rds"))),
       Exp = readRDS(file.path(work_dir, paste0("Exp_clu_",     project_name, ".rds"))))
}
cond_Ctrl <- load_condition(proj[["Control"]])
cond_AD   <- load_condition(proj[["AD"]])

tryCatch({
  ## (i) cell type -> cell type interaction strength, AD minus Control
  W_c <- cond_Ctrl$LR@net$weight
  W_a <- cond_AD$LR@net$weight
  ct  <- intersect(rownames(W_c), rownames(W_a))
  long <- as.data.frame(as.table(W_a[ct, ct] - W_c[ct, ct]), stringsAsFactors = FALSE)
  colnames(long) <- c("Sender", "Receiver", "AD_minus_Control")
  long$Control <- as.vector(W_c[ct, ct])
  long$AD      <- as.vector(W_a[ct, ct])
  long <- long[order(-abs(long$AD_minus_Control)), ]
  cat("\nLargest cell-type-to-cell-type interaction-strength differences (AD - Control):\n")
  print(head(long, 10))
  write.csv(long, file.path(table_dir, "CCI_celltype_weights_AD_vs_Control.csv"), row.names = FALSE)

  ## (ii) ligand -> receptor pairs signaling INTO the receiver cell type
  lr_c <- as.data.frame(Extract_LR_Prob(cond_Ctrl$LR, target_type = target_type))
  lr_a <- as.data.frame(Extract_LR_Prob(cond_AD$LR,   target_type = target_type))
  if (all(c("From", "To", "Weight") %in% names(lr_c)) && all(c("From", "To", "Weight") %in% names(lr_a))) {
    agg <- function(x) x %>% dplyr::group_by(From, To) %>%
      dplyr::summarise(Weight = sum(Weight), .groups = "drop")
    lr_cmp <- dplyr::full_join(agg(lr_c), agg(lr_a), by = c("From", "To"),
                               suffix = c("_Control", "_AD"))
    lr_cmp$Weight_Control[is.na(lr_cmp$Weight_Control)] <- 0
    lr_cmp$Weight_AD[is.na(lr_cmp$Weight_AD)]           <- 0
    lr_cmp$AD_minus_Control <- lr_cmp$Weight_AD - lr_cmp$Weight_Control
    lr_cmp <- lr_cmp[order(-abs(lr_cmp$AD_minus_Control)), ]
    cat("\nLigand->Receptor pairs signaling into ", target_type,
        " with the largest weight difference (AD - Control):\n", sep = "")
    print(head(as.data.frame(lr_cmp), 15))
    write.csv(lr_cmp, file.path(table_dir, paste0("LR_pairs_into_", target_type, "_AD_vs_Control.csv")), row.names = FALSE)
  } else {
    cat("\nExtract_LR_Prob output did not have From/To/Weight columns - skipping the L-R table.\n")
  }
}, error = function(e) {
  cat("\nSection 10b (numeric CellChat comparison) failed: ", conditionMessage(e),
      "\nThis does not affect anything else - continuing to the figures.\n", sep = "")
})

## --- 11. Helper needed by Section 12 ----------------------------------------
## A target gene's Section-6 DE status is not enough to plot it: PlotXT_HeatMap
## pivots pathways into a Receptor x SSC matrix that collapses if a gene's
## pathways all share a single Receptor or a single SSC. Walk down the ranked
## list and use the first gene with at least 2 distinct values in BOTH.
find_plottable_gene <- function(results, ranked_genes, min_unique = 2) {
  for (g in ranked_genes) {
    sub <- dplyr::filter(results, Target == g)
    if (dplyr::n_distinct(sub$Receptor) >= min_unique &&
        dplyr::n_distinct(sub$SSC) >= min_unique) {
      return(g)
    }
  }
  NA_character_
}

## ============================================================================
## --- 12. The full SigXTalk visualization gallery, for BOTH conditions -----
## Source: https://github.com/LithiumHou/SigXTalk/blob/master/vignettes/demo.md
## ("Visualize the results"), using this pipeline's own object names.
## Every figure is wrapped so that one failing figure is reported and skipped
## instead of stopping the rest of the gallery or the other condition.
## ============================================================================

try_step <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    cat("\n[skipped] ", label, " failed: ", conditionMessage(e), "\n", sep = "")
    NULL
  })
}

## Open a PNG, evaluate the (lazy) plotting expression inside it, always close.
save_base_png <- function(path, expr, width = 1600, height = 1600, res = 200) {
  png(path, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  invisible(NULL)
}

run_section12_viz <- function(condition, LR_original, Exp_clu, results_now, target_type, work_dir) {

  fig_dir <- file.path(work_dir, "figures", paste0("GSE174367_", condition, "_full"))
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\n=== Section 12 for", condition, "- figures will be written to:", fig_dir, "===\n")

  ## --- 12a. Visualization of cell-cell communication -----------------------
  ## (i) Overview: aggregated interaction strength between every pair of cell types.
  ## NOTE: circle colors follow each object's own factor order, so compare
  ## conditions by the text labels, never by color. Edge widths are scaled within
  ## each plot - use the tables from Section 10b for real magnitudes.
  try_step("12a CCI overview circle", {
    groupSize <- as.numeric(table(LR_original@idents))
    save_base_png(file.path(fig_dir, "12a_CCI_overview_circle.png"), {
      CellChat::netVisual_circle(LR_original@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                                  label.edge = FALSE, title.name = "Interaction weights/strength",
                                  vertex.label.cex = 2)
    })
  })

  ## (ii) Which ligand-receptor pairs from OTHER cell types signal INTO the receiver?
  try_step("12a L-R circle into receiver", {
    CCC_threshold <- 0.1
    LR_Pairprob <- Extract_LR_Prob(LR_original, target_type = target_type)
    LR_Pairprob <- LR_Pairprob[LR_Pairprob$Weight >= CCC_threshold * max(LR_Pairprob$Weight), ]
    cat("\nL-R pairs signaling into", target_type, "above threshold:", nrow(LR_Pairprob), "\n")
    save_base_png(file.path(fig_dir, paste0("12a_LR_circle_", target_type, ".png")), {
      PlotCCI_CirclePlot(LR_Pairprob, topk = 30)
    })
  })

  ## --- 12b. Overview of the crosstalk patterns -----------------------------
  ## $outer: histogram of how many crosstalk pathways explain each target gene
  ## $inner: bar chart of the top genes ranked by that count
  try_step("12b PlotXT_Counts", {
    ps <- PlotXT_Counts(results_now, top_percent = 1)
    ggsave(file.path(fig_dir, "12b_crosstalk_counts_distribution.png"), ps$outer, width = 7, height = 5)
    ggsave(file.path(fig_dir, "12b_crosstalk_counts_topgenes.png"),     ps$inner, width = 7, height = 6)
  })

  ## --- 12c. Genes with the most crosstalk pathways -------------------------
  CC_pair_results <- Aggregate_Causality(results_now, data_type = "Target")
  Counts_pathway  <- Count_Crosstalk(results_now, KeyGenes = NULL, data_type = "Target", verbose = FALSE)
  TG_ranked       <- sort(Counts_pathway, decreasing = TRUE) %>% names()
  TG_top15        <- TG_ranked[1:min(15, length(TG_ranked))]
  cat("\nTop target genes by number of crosstalk pathways (", condition, "):\n", sep = "")
  print(head(Counts_pathway[TG_ranked], 15))
  write.csv(data.frame(Target = TG_ranked, n_pathways = as.numeric(Counts_pathway[TG_ranked])),
            file.path(work_dir, "tables", paste0("pathway_counts_per_target_", condition, ".csv")),
            row.names = FALSE)

  p_rectg <- tryCatch(
    PlotXT_RecTGHeatmap(CC_pair_results, Exp_clu, KeyTG = TG_top15, topk = 100),
    error = function(e) {
      cat("\nPlotXT_RecTGHeatmap failed (", conditionMessage(e), ") - the top pathways for\n",
          "these ", length(TG_top15), " genes collapsed to a single Receptor or Target,\n",
          "which this function's internal reshape doesn't guard against.\n",
          "Skipping this figure for ", condition, ".\n", sep = "")
      NULL
    }
  )
  if (!is.null(p_rectg)) {
    ggsave(file.path(fig_dir, "12c_receptor_target_heatmap.png"), p_rectg, width = 10, height = 10)
  }

  ## --- 12d. Detailed crosstalk pattern of specific genes -------------------
  gene_for_heatmap <- find_plottable_gene(results_now, TG_ranked)

  if (!is.na(gene_for_heatmap)) {
    try_step("12d PlotXT_HeatMap", {
      p_heat <- PlotXT_HeatMap(results_now, gene_for_heatmap, "TG")   # "TG" == "Target"
      save_base_png(file.path(fig_dir, paste0("12d_heatmap_", gene_for_heatmap, ".png")), {
        ComplexHeatmap::draw(p_heat)
      }, width = 2200, height = 2200)
    })
  } else {
    cat("\nNo target gene has >=2 distinct Receptors/SSCs in this condition -",
        "skipping PlotXT_HeatMap (a real sparsity finding, not a bug).\n")
    gene_for_heatmap <- TG_top15[1]
  }

  ## Alluvial (river) plot of Receptor -> SSC -> Target flow for the top-3 genes
  alluvial_genes <- TG_top15[1:min(3, length(TG_top15))]
  p_alluvial <- tryCatch(
    PlotXT_Alluvial(results_now, alluvial_genes, min_weight = 0.8),
    error = function(e) {
      cat("\nPlotXT_Alluvial with min_weight = 0.8 failed (", conditionMessage(e),
          ") - retrying with min_weight = 0 (no weight filter).\n", sep = "")
      tryCatch(PlotXT_Alluvial(results_now, alluvial_genes, min_weight = 0),
               error = function(e2) {
                 cat("[skipped] 12d alluvial: ", conditionMessage(e2), "\n", sep = "")
                 NULL
               })
    }
  )
  if (!is.null(p_alluvial)) {
    ggsave(file.path(fig_dir, "12d_alluvial_top3.png"), p_alluvial, width = 8, height = 6)
  }

  ## Fidelity (how selectively each Receptor/SSC pair regulates THIS gene) and
  ## Specificity (how selectively vs. every other target)
  try_step("12d PlotXT_FidSpe", {
    p_fidspe <- PlotXT_FidSpe(results_now, gene_for_heatmap, threshold = 0.01)
    ggsave(file.path(fig_dir, paste0("12d_fidelity_specificity_", gene_for_heatmap, ".png")),
           p_fidspe, width = 10, height = 6)
  })

  cat("\nSection 12 finished for", condition, "- figures are in:", fig_dir, "\n")
  invisible(NULL)
}

## --- Run it for Control, then AD - each from its OWN checkpointed data ----
run_section12_viz("Control", cond_Ctrl$LR, cond_Ctrl$Exp, results_Ctrl, target_type, work_dir)
run_section12_viz("AD",      cond_AD$LR,   cond_AD$Exp,   results_AD,   target_type, work_dir)

## Record the package versions this run used
writeLines(capture.output(sessionInfo()), file.path(work_dir, "sessionInfo_full_run.txt"))

cat("\n\nDone (all cells, no downsampling).\n",
    "Figures: ", file.path(work_dir, "figures"), "\n",
    "Tables:  ", table_dir, "\n", sep = "")
