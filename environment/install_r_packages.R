## Installs the R packages used by analysis/01_sigxtalk_pipeline.Rmd and analysis/02_de_direction_check.Rmd.
## Run once, in R >= 4.1 (the reported run used R 4.6.1; see results/sessionInfo_full_run.txt for versions).

cran_packages <- c(
  "devtools", "BiocManager", "rmarkdown", "knitr",
  "Seurat", "dplyr", "tidyr", "ggplot2", "patchwork", "ggalluvial",
  "foreach", "doSNOW", "snow", "ranger", "tidymodels", "igraph", "hdf5r"
)
install.packages(setdiff(cran_packages, rownames(installed.packages())))

## Bioconductor dependencies of SigXTalk / CellChat
BiocManager::install(c("Biobase", "BiocNeighbors", "ComplexHeatmap", "BiocGenerics"))

## CellChat (development version from GitHub; the reported run used 2.2.0.9001)
devtools::install_github("jinworks/CellChat")

## SigXTalk R package (installs as "SigXTalkR"); instructions from https://github.com/LithiumHou/SigXTalk
devtools::install_github("LithiumHou/SigXTalk", dependencies = TRUE, upgrade = "always")
