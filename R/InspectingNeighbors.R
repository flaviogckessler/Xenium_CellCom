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
###############################
# Working with polygon geometry
nodes <- gs_crop@nodes
nodes <- nodes %>% select(c(vertex,x,y,nCount_Xenium,nFeature_Xenium,geometry))

st_geometry_type(nodes$geometry) %>% table()

###############################
ggplot(data = st_as_sf(nodes[1:9,]))+
  geom_sf()+
  #geom_point(aes(x=x,y=y))+
  geom_text(aes(x=x,y=y,label = vertex))+
  theme_minimal()

intersect_list <- st_intersects(st_as_sf(nodes[1:9,]),remove_self=TRUE)
intersect_list

neighbors_df <- st_join(st_as_sf(nodes[1:9,]), st_as_sf(nodes[1:9,]), join = st_intersects,remove_self=TRUE)
neighbors_df$status <- ifelse(is.na(neighbors_df$vertex.y), "No", "Yes")

ggplot(data = neighbors_df)+
  geom_sf(aes(fill=status))+
  geom_text(aes(x=x.x,y=y.x,label = vertex.x))+
  scale_fill_manual(values = c("Yes"="#2b8cbe","No"="grey"),
                    na.value = "white")+
  theme_minimal()+
  labs(fill="Have Neighbor?",x="x",y="y")

buffered_nodes <- st_buffer(st_as_sf(nodes[1:9,]),dist = 0.0005)
ggplot(data = buffered_nodes)+
  geom_sf()+
  #geom_point(aes(x=x,y=y))+
  geom_text(aes(x=x,y=y,label = vertex))+
  theme_minimal()

neighbors_df <- st_join(buffered_nodes, buffered_nodes, join = st_intersects,remove_self=TRUE)
st_intersects(buffered_nodes,remove_self=TRUE)
neighbors_df$status <- ifelse(is.na(neighbors_df$vertex.y), "No", "Yes")

ggplot(data = neighbors_df)+
  geom_sf(aes(fill=status))+
  geom_text(aes(x=x.x,y=y.x,label = vertex.x))+
  scale_fill_manual(values = c("Yes"="#2b8cbe","No"="grey"),
                    na.value = "white")+
  theme_minimal()+
  labs(fill="Have Neighbor?",x="x",y="y")
##########
nodes <- st_as_sf(nodes[1:9,])
neighbor_list <- st_is_within_distance(nodes,nodes, dist = 0.0005,remove_self=TRUE)
neighbor_list
length(neighbor_list)
neighbor_list[[1]]
length(neighbor_list[[1]])

neighbor_list <- as.data.frame(neighbor_list)
View(neighbor_list)
##########
# Assigning spatial metrics as nodes features
gs_vertex_attr(gs_crop,"geometry")
# define function that calculate the max distance polygon diameter
maxDistPol <- function(polygon){
  # Extract polygon points
  polygon_points <- st_cast(polygon, "POINT")
  polygon_points <- polygon_points[-1]
  
  # Find the absolute maximum distance between any two vertices
  #max_distance <- max(st_distance(polygon_points, polygon_points[-1]))
  max_distance <- max(st_distance(polygon_points, polygon_points))
  
  return(max_distance)
}

# Assign maximum distance (polygon diameter)
all_geom <- gs_crop@nodes %>% select(vertex,geometry)

for(i in 1:length(nodes$vertex)){
  all_geom$maxDist[i] <- maxDistPol(all_geom$geometry[i])
}