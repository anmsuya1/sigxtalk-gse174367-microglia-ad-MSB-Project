# Data

## Source

Single-nucleus RNA-seq of human prefrontal cortex from Alzheimer's disease (AD) and control donors.

* Study: Morabito S. et al., *Single-nucleus chromatin accessibility and transcriptomic characterization of
  Alzheimer's disease*, Nature Genetics 53, 1143-1155 (2021)
* GEO accession: **GSE174367** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367>

## Files the pipeline expects in this folder

| File | Size | In this repository? | MD5 checksum |
|---|---:|---|---|
| `GSE174367_snRNA-seq_cell_meta.csv` | 4.5 MB | yes | `af4f2c63a1c6683efce536cd6698e252` |
| `GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5` | 274 MB | **no** (exceeds GitHub's 100 MB file limit) | `ef6ca02634989e2933f010ea6d6926c5` |

## How to get the count matrix

1. Open the GEO page above and go to the *Supplementary file* table at the bottom.
2. Download `GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5` (the snRNA-seq filtered feature-barcode matrix).
   If GEO serves it compressed (`.h5.gz`), decompress it so that the file ends in `.h5`.
3. Place it in this folder (`data/`), next to the metadata file.
4. Optional check: `md5sum GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5` (macOS: `md5 <file>`) should match
   the checksum in the table.

The metadata file is also available from the same GEO page (you can replace the included copy; the checksum above is for
the copy used in the reported run).

## Metadata columns

`Barcode`, `SampleID`, `Diagnosis` (AD / Control), `Batch`, `Cell.Type` (ASC, EX, INH, MG, ODC, OPC, PER.END),
`cluster`, `Age`, `Sex`, `PMI`, `Tangle.Stage`, `Plaque.Stage`, `RIN`.

The count matrix has 61,770 barcodes; the metadata annotates 61,472 of them. The pipeline keeps the annotated nuclei.
