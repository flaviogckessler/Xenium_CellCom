# https://htmlpreview.github.io/?https://github.com/jinworks/SpatialCellChat/blob/master/tutorial/SpatialCellChat_analysis_of_spatial_transcriptomics_data.html

# Install helper packages if needed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("SpatialCellChat", quietly = TRUE)) {
  devtools::install_github("jinworks/SpatialCellChat")
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  devtools::install_github("jinworks/CellChat")
}

if (!requireNamespace("BiocNeighbors", quietly = TRUE)) {
  BiocManager::install("BiocNeighbors")
}

if (!requireNamespace("MERINGUE", quietly = TRUE)) {
  devtools::install_github("JEFworks-Lab/MERINGUE")
}

if (!requireNamespace("ALRA", quietly = TRUE)) {
  devtools::install_github("KlugerLab/ALRA")
}

if (!requireNamespace("RcppML", quietly = TRUE)) {
  devtools::install_github("zdebruine/RcppML")
}

if (!requireNamespace("presto", quietly = TRUE)) {
  devtools::install_github("immunogenomics/presto")
}

# Load packages
library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratObject)
library(SeuratData)
library(dplyr)
library(tibble)
