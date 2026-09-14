# reference tutorial: https://sysbiolab.github.io/RGraphSpace/articles/nested-geometries.html

# Install packages
if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
if (!require("SpatialFeatureExperiment", quietly = TRUE)){BiocManager::install("SpatialFeatureExperiment")}
if (!require("RGraphSpace", quietly = TRUE)){
  if (!require("remotes", quietly = TRUE)){install.packages("remotes")}
  remotes::install_github("sysbiolab/RGraphSpace", build_vignettes=TRUE)}
if (!require("sf", quietly = TRUE)){install.packages("sf")}
if (!require("terra", quietly = TRUE)){install.packages("terra")}

# Check versions
if (packageVersion("RGraphSpace") < "1.5.4"){
  message("Need to update 'RGraphSpace' for this vignette")
  remotes::install_github("sysbiolab/RGraphSpace")}

# Load packages
library("RGraphSpace")
library("SpatialFeatureExperiment")
library("sf")
library("terra")
library("patchwork")
#################
# Batch:
# Download output files in a 'localdir' directory
#wget https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_Human_Colon_Cancer_P2_CRC_Add_on_FFPE/Xenium_V1_Human_Colon_Cancer_P2_CRC_Add_on_FFPE_outs.zip

# Extract the outputs
#unzip Xenium_V1_Human_Colon_Cancer_P2_CRC_Add_on_FFPE_outs.zip
#################