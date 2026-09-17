library("RGraphSpace")
library("sf")
library("ggplot2")
library("deldir")
library("igraph")

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
    xlab = "Tissue coordinate 1",
    ylab = "Tissue coordinate 2")

p2
###############################
colnames(gs_crop@nodes)
gs_vertex_attr(gs_crop,"x")[1:5]
gs_crop@nodes$x[1:5]
gs_vertex_attr(gs_crop,"y")[1:5]
gs_crop@nodes$y[1:5]
gs_vertex_attr(gs_crop,"nodeLabel")[1:5]


triangulation <- deldir(gs_crop@nodes$x,gs_crop@nodes$y)
# Extract the edges (connections between points)
edges_matrix <- triangulation$delsgs[, c("ind1", "ind2")]
head(edges_matrix)

edges_matrix$ind1 <- gs_vertex_attr(gs_crop,"name")[edges_matrix$ind1]
edges_matrix$ind2 <- gs_vertex_attr(gs_crop,"name")[edges_matrix$ind2]
colnames(edges_matrix) <- c("from","to")
head(edges_matrix)

# Adding edges to the GraphSpace object
gs_crop_dnay <- gs_crop |> gs_add_edges(edges_matrix)

# Teste if the number of edges in the edge_matrix is equal to the number of edges in the GraphSpace object
length(edges_matrix$from) == length(gs_crop_dnay@edges$vertex1)

gs_edge_attr(gs_crop_dnay, "arrowType") <- 0

p3 <- ggplot(gs_crop_dnay) +
  geom_edgespace()+
  geom_nodespace(
    mapping = aes(colour = log2(Slc17a7 + 1)),
    size = .5,pch = 19) +
  scale_colour_continuous(palette = cpal, limits = data_range) +
  theme_gspace_coords(
    theme = "th3",is_norm = TRUE,
    xlab = "Coordenada do tecido 1",
    ylab = "Coordenada do tecido 2")

p3

p1 <- ggplot(gs_crop) +
  geom_sf(mapping = aes(geometry = geometry)) +
  scale_fill_continuous(palette = cpal, limits = data_range) +
  geom_nodespace(colour = "black", size = 0.2, pch = 19) +
  theme_gspace_coords(theme = "th3",is_norm = TRUE)+
  theme_minimal()+
  labs(title = "Cells segmentation",x=NULL,y=NULL)

p2 <- ggplot(gs_crop_dnay)+
  geom_nodespace(color="black",
                 size = .5,pch = 19) +
  scale_colour_continuous(palette = cpal, limits = data_range) +
  theme_gspace_coords(theme = "th3",is_norm = TRUE)+
  theme_minimal()+
  labs(title = "Cells centroids",x=NULL,y=NULL)

p3 <- ggplot(gs_crop_dnay) +
  geom_edgespace()+
  geom_nodespace(color="black",
    size = .5,pch = 19) +
  scale_colour_continuous(palette = cpal, limits = data_range)+
  theme_gspace_coords(theme = "th3",is_norm = TRUE)+
  theme_minimal()+
  labs(title = "Delaunay triangulation",x=NULL,y=NULL)

p1 + p2 + p3

ggsave("images/GrapsSpace_Delaunay.png",p1 + p2 + p3,
       height = 5,width = 12,units = "in",dpi = 120)
######################################################
# Calculate the greates distance btween any pair of its vertices
# load function that calculate the maximum distance of a polygon (polygon diameter)
source("R/maxDistPol_func.R")

# Assign maximum distance (polygon diameter)
all_geom <- gs_crop@nodes %>% select(vertex,geometry)

for(i in 1:length(all_geom$geometry)){
  all_geom$maxDist[i] <- maxDistPol(all_geom$geometry[i])
}
gs_vertex_attr(gs_crop_dnay, "maxDist") <- all_geom$maxDist
rm(all_geom)

# Store the meanMaxDist between all the cell polygons
meanMaxDist <- mean(gs_vertex_attr(gs_crop_dnay, "maxDist"))
######################################################
# Finding neighbors by cell segmentation
neighbor_list <- st_is_within_distance(st_as_sf(gs_crop_dnay@nodes),
                                       st_as_sf(gs_crop_dnay@nodes),
                                       dist = 0.0005,
                                       remove_self=TRUE)
neighbor_list
####
neighbor_df<- as.data.frame(neighbor_list)
head(neighbor_df)

neighbor_df$row.id <- gs_vertex_attr(gs_crop,"name")[neighbor_df$row.id]
neighbor_df$col.id <- gs_vertex_attr(gs_crop,"name")[neighbor_df$col.id]
colnames(neighbor_df) <- c("from","to")
head(neighbor_df)

# Adding edges to the GraphSpace object by cell-cell proximity
gs_crop_prox <- gs_crop |> gs_add_edges(neighbor_df)

p4 <- ggplot(gs_crop_prox) +
  geom_edgespace()+
  geom_nodespace(color="black",
                 size = .5,pch = 19) +
  scale_colour_continuous(palette = cpal, limits = data_range)+
  theme_gspace_coords(theme = "th3",is_norm = TRUE)+
  theme_minimal()+
  labs(title = "Cell proximity",x=NULL,y=NULL)
p4

gs_crop_dnay@graph
gs_crop_prox@graph

g <- ggpubr::ggarrange(p1,p2,p3,p4,
                       ncol = 2, nrow = 2,
                       labels = c("A", "B", "C", "D"))
g

ggsave("images/GrapsSpace_Delaunay_Proximity.png",g,
       height = 10,width = 10,units = "in",dpi = 120)
