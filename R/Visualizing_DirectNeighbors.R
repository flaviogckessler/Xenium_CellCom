library("RGraphSpace")
library("sf")
library("ggplot2")
library("igraph")
library("dplyr")

#install.packages("gifski")
library(gifski)

# download crop
diretorio_dados <- "data"
url_dados <- paste0(
  "https://raw.githubusercontent.com/flaviogckessler/Curso_ST/",
  "main/data/Xenium_GS_Crop.RDS")

arquivo_dados <- file.path(diretorio_dados, "Xenium_GS_Crop.RDS")

if (!file.exists(arquivo_dados)) {
  utils::download.file(url_dados,arquivo_dados,mode = "wb")}

gs_crop <- readRDS(arquivo_dados)

if (!inherits(gs_crop, "GraphSpace")) {
  stop("O arquivo baixado não contém um objeto da classe GraphSpace.")}

gs_crop
#########
# Set color palette and data range for use across plots
cpal <- hcl.colors(100, palette = "Spectral", rev = T)
data_range <- c(0, 7)
# Set a reusable theme
my_theme <- theme_gspace_coords(theme = "th3", is_norm = TRUE, 
                                xlab = "Tissue coordinates 1", ylab = "Tissue coordinates 2")
#########
# load function that calculate that find the spot/cells neighbors based on the cell segmentation info
source("R/gs_neighbors_func.R")
#########
p1 <- ggplot(gs_crop) +
  geom_sf(
    mapping = aes(
      geometry = geometry,
      fill = log2(Slc17a7 + 1))) +
  scale_fill_continuous(palette = cpal, limits = data_range) +
  geom_nodespace(colour = "black", size = 0.2, pch = 19) +
  my_theme

p1
###############################
gs_neig <- gs_neighbors(gs_crop,dist=0.0005)

p2 <- ggplot(gs_neig) +
  geom_edgespace(color="red")+
  geom_nodespace(color="black",
                 size = .2,pch = 19) +
  scale_colour_continuous(palette = cpal, limits = data_range)+
  theme_gspace_coords(theme = "th3",is_norm = TRUE)+
  my_theme

p2

p3 <- ggplot(gs_neig) +
  geom_sf(mapping = aes(geometry = geometry,fill = log2(Slc17a7 + 1))) +
  scale_fill_continuous(palette = cpal, limits = data_range) +
  #geom_nodespace(colour = "black", size = 0.2, pch = 19) +
  geom_edgespace(color="red",alpha=0.3)+
  my_theme

p3
#####################
# Crop to the smaller region and re-normalize coordinates
gs_neig2 <- cropGraphSpace(gs_neig, 
                           xmin = 0.45, xmax = 0.6, 
                           ymin = 0.2, ymax = 0.35)
gs_neig2 <- normalizeGraphSpace(gs_neig2, mar = 0, norm.geometry = TRUE)

ggplot(gs_neig2) +
  geom_sf(mapping = aes(geometry = geometry,fill = log2(Slc17a7 + 1))) +
  scale_fill_continuous(palette = cpal, limits = data_range) +
  #geom_nodespace(colour = "black", size = 0.2, pch = 19) +
  geom_edgespace(color="red",alpha=0.8)+
  my_theme


gs_edge_attr(gs_neig2,"arrowType") <- "13|-|13"

ggplot(gs_neig2) +
  geom_sf(mapping = aes(geometry = geometry,fill = log2(Slc17a7 + 1))) +
  scale_fill_continuous(palette = cpal, limits = data_range) +
  geom_edgespace(color="black",alpha=0.8)+
  my_theme

ggplot(gs_neig2) +
  geom_sf(mapping = aes(geometry = geometry,fill = log2(Slc17a7 + 1))) +
  scale_fill_continuous(palette = cpal, limits = data_range) +
  geom_nodespace(colour = "black", size = 15, pch = 19,alpha=0) +
  geom_edgespace(color="black",alpha=0.8)+
  my_theme

library("patchwork")
glyph_list()
