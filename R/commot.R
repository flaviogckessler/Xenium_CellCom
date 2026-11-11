# Load packages
library("RGraphSpace")
library("SpatialExperiment")
library("patchwork")
library("sf")

#BiocManager::install("OSTA.data")
library("OSTA.data")
#BiocManager::install("VisiumIO")
library("VisiumIO")
################################################################################
### Loading the dataset
################################################################################
# Set a temporary directory for download
id <- "VisiumHD_HumanColon_Oliveira"
bfc <- BiocFileCache::BiocFileCache()
qid <- BiocFileCache::bfcquery(bfc, id)$rid
BiocFileCache::bfcremove(bfc, qid)

dir.create(td <- tempfile())

# Retrieve data from OSF repo
unzip(OSTA.data_load(id), exdir = td)

# Read into 'SpatialExperiment'
spe <- TENxVisiumHD(
  segmented_outputs = file.path(td, "segmented_outputs"),
  format = "h5", 
  images = "hires") |>
  import()

# Assign unique symbols to rownames
rownames(spe) <-  make.unique(rowData(spe)$Symbol)

library(scater)
plotExpression(
  spe, 
  features = c("GAPDH", "ACTB"), 
  x = "label", 
  assay.type = "logcounts"
)
################################################################################
### Setup Python from R
################################################################################

#-------------------------------------------------------------------------------
# Option 1, a named virtual environment (any reticulate version)
# Restart R before running, so Python isn't already loaded
#-------------------------------------------------------------------------------

# Load reticulate, the R <-> Python bridge
library(reticulate)

# Install Python 3.10 (skipped if already installed)
py_path <- install_python(version = "3.10:latest")

# Environment name
env <- "commot_env"

# Start clean: remove any previous environment with this name
if (virtualenv_exists(env)) virtualenv_remove(env, confirm = FALSE)

# Create the environment with NumPy 1.x (COMMOT breaks on NumPy 2)
virtualenv_create(envname = env, python = py_path, packages = "numpy<2")

# Install numpy 1.x; COMMOT uses features removed in NumPy 2.0
virtualenv_install(env, packages = "numpy<2")

# Install COMMOT, keeping NumPy 1.x and using a prebuilt POT
virtualenv_install(env, packages = c("numpy<2", "pot", "commot"),
  pip_options = "--only-binary=pot")

# Use the environment in this session
use_virtualenv(env, required = TRUE)

#-------------------------------------------------------------------------------
# Option 2, an automatic environment managed by reticulate;
# Requires reticulate >= 1.41;
# Restart R before running, so Python isn't already loaded
#-------------------------------------------------------------------------------

# Load reticulate, the R <-> Python bridge
library(reticulate)

# Declare Python and packages: NumPy 1.x (COMMOT breaks on NumPy 2),
# and no package releases after this date, so installs stay reproducible
# Declare COMMOT and NumPy 1.x (COMMOT breaks on NumPy 2), on the tested Python version
py_require(c("commot", "numpy<2"), python_version = "3.10")


################################################################################
### Run COMMOT
################################################################################

# Check Python is in use
py_config()

# Check that commot is importable (should print TRUE)
print(py_module_available("commot"))

# Import COMMOT into R (used to retrieve LR interactions from CellChatDB)
ct <- reticulate::import("commot")

# Retrieve LR interaction from CellChatDB
db <- ct$pp$ligand_receptor_database(
  species = "human", 
  database = "CellChat", 
  signaling_type = NULL)

names(db) <- c("ligand", "receptor", "pathway", "type")
head(db)

# interaction types &
# number of interactions
table(db$type)

# number of unique features
lr <- c(db$ligand, db$receptor)
ss <- strsplit(lr, "-")
gs <- sapply(ss, .subset, 1)
length(unique(unlist(gs)))

# filter for those fully covered by the panel
keep <- apply(db, 1, \(.) {
  rs <- strsplit(.["receptor"], "_")
  lr <- c(.["ligand"], unlist(rs))
  all(lr %in% rownames(spe))
})
(db <- db[keep, ])

# subset data to features being considered
ss <- strsplit(db$receptor, "_")
rs <- sapply(ss, .subset, 1)
sub <- spe[unique(c(db$ligand, rs)), ]

# assign spatial coordinates to 'spatial' 
# reduced dimension slot (required by 'COMMOT')
xy <- spatialCoords(sub)
reducedDim(sub, "spatial") <- xy

# get indices of cells that fall
# into a (2mm) x (2mm) window
cs <-
  xy[, 1] >= 1500 & xy[, 1] <= 3500 &
  xy[, 2] >= 2000 & xy[, 2] <= 4000

# print number of genes x cells left
dim(sub <- sub[, cs]) 

# realize expression matrix
library("Matrix")
logcounts(sub) <- as(logcounts(sub), "dgCMatrix")

# convert 'SpatialExperiment' (R) to 
# 'AnnData' (Python) using 'anndataR'
ad <- as_AnnData(sub, 
                 x_mapping="logcounts", 
                 output_class="ReticulateAnnData")

# run CCC inference using 'COMMOT'
ct <- import("commot")
ct$tl$spatial_communication(ad,
                            database_name="CellChatDB",
                            dis_thr=10,  
                            df_ligrec=db,
                            heteromeric=TRUE,
                            pathway_sum=TRUE,
                            heteromeric_rule="min",
                            heteromeric_delimiter="_")

# extract 'data.frame's containing
# sender & receiver CCC estimates
ccc <- list(
  s=ad$obsm$`commot-CellChatDB-sum-sender`,
  r=ad$obsm$`commot-CellChatDB-sum-receiver`)

# add CCC inference results as cell metadata
for (df in ccc) colData(sub)[names(df)] <- df