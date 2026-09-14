library("RGraphSpace")
library("Seurat")
library("SeuratObject")
library("sf")
library("ggplot2")
library("patchwork")
library("dplyr")

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
# Inspect the data range
# log2(range(gs[["fdata"]]) + 1)

# Set color palette and data range for use across plots
cpal <- hcl.colors(100, palette = "Spectral", rev = T)
data_range <- c(0, 7)
# Set a reusable theme
my_theme <- theme_gspace_coords(theme = "th3", is_norm = TRUE, 
                                xlab = "Tissue coordinates 1", ylab = "Tissue coordinates 2")
#########
p1 <- ggplot(gs_crop) +
  geom_nodespace(
    mapping = aes(colour = log2(Slc17a7 + 1)),
    size = 1.5,pch = 19) +
  scale_colour_continuous(palette = cpal, limits = data_range) +
  theme_gspace_coords(
    theme = "th3",is_norm = TRUE,
    xlab = "Coordenada do tecido 1",
    ylab = "Coordenada do tecido 2")

p1
p2 <- ggplot(gs_crop) +
  geom_sf(
    mapping = aes(
      geometry = geometry,
      fill = log2(Slc17a7 + 1))) +
  scale_fill_continuous(palette = cpal, limits = data_range) +
  geom_nodespace(colour = "black", size = 0.2, pch = 19) +
  theme_gspace_coords(
    theme = "th3",is_norm = TRUE,
    xlab = "Coordenada do tecido 1",
    ylab = "Coordenada do tecido 2")

p2
###############################
# Working with polygon geometry
nodes <- gs_crop@nodes
nodes <- nodes %>% select(c(vertex,x,y,nCount_Xenium,nFeature_Xenium,geometry))

c <- ggplot(nodes, aes(x="",y=nCount_Xenium))+
  geom_violin()+
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5)+
  labs(x = NULL)
f <- ggplot(nodes, aes(x="",y=nFeature_Xenium))+
  geom_violin()+
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5)+
  labs(x = NULL)

st_geometry_type(nodes$geometry) %>% table()
####
# 1. Calculate the area

nodes$area <- st_area(nodes$geometry)
a <- ggplot(nodes, aes(x="",y=area))+
  geom_violin()+
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5)+
  labs(x = NULL)
####
# 2. Calculate the max distance polygon diameter
# Define function
maxDistPol <- function(polygon){
  # Extract polygon points
  polygon_points <- st_cast(polygon, "POINT")
  polygon_points <- polygon_points[-1]
  
  # Find the absolute maximum distance between any two vertices
  #max_distance <- max(st_distance(polygon_points, polygon_points[-1]))
  max_distance <- max(st_distance(polygon_points, polygon_points))
  
  return(max_distance)
}

for(i in 1:length(nodes$vertex)){
  nodes$maxDist[i] <- maxDistPol(nodes$geometry[i])
}

md <- ggplot(nodes, aes(x="",y=maxDist))+
  geom_violin()+
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5)+
  labs(x = NULL)
####
# 3. Calculate the perimeter
nodes$perimeter <- st_perimeter(nodes$geometry)

p <- ggplot(nodes, aes(x="",y=perimeter))+
  geom_violin()+
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5)+
  labs(x = NULL)
####
a + p + md + plot_annotation(title = "Spatial Metrics from cells")
#c + f + a + md + p + plot_annotation(title = "Spatial Metrics from cells")
##########
# Assigning spatial metrics as nodes features
gs_vertex_attr(gs_crop,"geometry")

class(gs_crop@nodes)
# Assgin area
gs_vertex_attr(gs_crop, "area") <- st_area(gs_crop@nodes$geometry)
# Assign maximum distance (polygon diameter)
all_geom <- gs_crop@nodes %>% select(vertex,geometry)

for(i in 1:length(nodes$vertex)){
  all_geom$maxDist[i] <- maxDistPol(all_geom$geometry[i])
}
gs_vertex_attr(gs_crop, "maxDist") <- all_geom$maxDist
rm(all_geom)
# Assign perimeter
gs_vertex_attr(gs_crop, "perimeter") <- st_perimeter(gs_crop@nodes$geometry)

a <- ggplot(gs_crop@nodes, aes(x="",y=area))+
  geom_violin()+
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5)+
  labs(x = NULL)

md <- ggplot(gs_crop@nodes, aes(x="",y=maxDist))+
  geom_violin()+
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5)+
  labs(x = NULL, y="Maximum Distance - Cell Diameter")

p <- ggplot(gs_crop@nodes, aes(x="",y=perimeter))+
  geom_violin()+
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5)+
  labs(x = NULL)
####
spatmetri <- a + p+ md + plot_annotation(title = "Spatial Metrics from cells")
ggsave("images/SpatialFeatures_RawCrop.png",spatmetri,
       height = 4,width = 5,units = "in",dpi = 120)
##############################################
# Gene expression correlation with spatial metrics
teste<-gs_crop@fdata
colSums(teste)
teste[,"Slc17a7"]

#BiocManager::install("microbiome")
#library(microbiome)
install.packages("psych")
library("psych")

# rows = nodes/spots
# columns = genes
expr <- as.data.frame(gs_crop@fdata)
# rows = same nodes/spots
# columns = numeric metadata features
meta <- gs_crop@nodes %>% select(c(area,maxDist,perimeter))

library(psych)
library(dplyr)
library(tidyr)
library(tibble)

# Pearson
pearson <- corr.test(
  x = expr,
  y = meta,
  method = "pearson",
  adjust = "BH",
  use = "pairwise"
)

# Spearman
spearman <- corr.test(
  x = expr,
  y = meta,
  method = "spearman",
  adjust = "BH",
  use = "pairwise"
)

pearson$r   # correlation coefficients
pearson$p   # adjusted P-values
pearson$n   # number of observations

spearman$r
spearman$p
spearman$n
###################
featureLevel <- c("area","perimeter","maxDist")
###################
result_pearson <- as.data.frame(pearson$r) |>
  rownames_to_column("gene") |>
  pivot_longer(
    -gene,
    names_to = "feature",
    values_to = "correlation") |>
  mutate(feature = factor(feature,levels = featureLevel))

pvalues_pearson <- as.data.frame(pearson$p) |>
  rownames_to_column("gene") |>
  pivot_longer(
    -gene,
    names_to = "feature",
    values_to = "p_adjusted") |>
      mutate(feature = factor(feature,levels = featureLevel))

result_pearson <- result_pearson |>
  left_join(
    pvalues_pearson,
    by = c("gene", "feature")) |>
  mutate(feature = factor(feature,levels = featureLevel))

######
result_spearman <- as.data.frame(spearman$r) |>
  rownames_to_column("gene") |>
  pivot_longer(
    -gene,
    names_to = "feature",
    values_to = "correlation"
  ) 

pvalues_spearman <- as.data.frame(spearman$p) |>
  rownames_to_column("gene") |>
  pivot_longer(
    -gene,
    names_to = "feature",
    values_to = "p_adjusted"
  )

result_spearman <- result_spearman |>
  left_join(
    pvalues_spearman,
    by = c("gene", "feature")
  )

###################
featureLevel <- c("area","perimeter","maxDist")

result_pearson <- result_pearson %>% 
  filter(p_adjusted < 0.01,abs(correlation)>=0.3)

result_pearson <- result_pearson %>% filter(abs(correlation) >= 0.3, p_adjusted < 0.05)

pearsonLevel <- result_pearson %>%
  select(gene,correlation)%>%
  arrange(correlation) %>%
  select(gene) %>% unique()

pearsonLevel <- pearsonLevel$gene

p <- microbiome::heat(as.data.frame(result_pearson), "feature", "gene", 
          fill = "correlation", limits = c(0,1),colours = c("white","red","darkred"),legend.text = "Pearson",
          star="p_adjusted",p.adj.threshold=0.01,association.threshold = 0.3) + 
  scale_x_discrete(limits = featureLevel)+
  scale_y_discrete(limits = pearsonLevel)+
  labs(title = "Gene expresion vs. spatial metrics",caption = "r >= 0.3") 

p 

ggsave("images/PearsonCorr_Gene_SpatialFeatures_RawCrop_Sig.png",p,
       height = 8,width = 5,units = "in",dpi = 120)

result_spearman <- result_spearman %>% filter(abs(correlation) >= 0.3, p_adjusted < 0.05)
spearmanLevel <- result_pearson %>%
  select(gene,correlation)%>%
  arrange(correlation) %>%
  select(gene) %>% unique()

spearmanLevel <- spearmanLevel$gene

g <- microbiome::heat(as.data.frame(result_spearman), "feature", "gene", 
                      fill = "correlation", limits = c(0,1),colours = c("white","red","darkred"),
                      legend.text = "Spearman",
                      star="p_adjusted",p.adj.threshold=0.01,association.threshold = 0.3) + 
  scale_x_discrete(limits = featureLevel)+
  scale_y_discrete(limits = spearmanLevel)+
  labs(title = "Gene expresion vs. spatial metrics",caption = "r >= 0.3") 

g 

ggsave("images/SpearmanCorr_Gene_SpatialFeatures_RawCrop_Sig.png",g,
       height = 8,width = 5,units = "in",dpi = 120)

#PRECISA
# 1. NÃO USAR CROP
# 2. NORMALIZAR DADOS PARA FAZER CORRELAÇÃO
##############################################
# Visual center and centroids
####

#point of view of polygons (Like point of inaccessbility)

# 1. Mathematical center (can fall outside the C-shape)
math_center <- st_centroid(my_polygon)

# 2. Visual center (guaranteed to be inside the C-shape)
visual_center <- st_point_on_surface(my_polygon)
