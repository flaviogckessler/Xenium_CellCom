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

# Inspecting the 'sf' geometry types
st_geometry_type(nodes$geometry) %>% table()
###############################
# Lets find the cell polygons which intersect with each other
# Viewing the cell polygons of the first nine cells:
ggplot(data = st_as_sf(nodes[1:9,]))+
  geom_sf()+
  geom_text(aes(x=x,y=y,label = vertex))+
  theme_minimal()
# If you get a closer look, you will see that even neighbor cells do not cross the polygons of all their neigbors (this will be handled dowstream)

# The simplest way to find possible neighbors is using the 'st-intersects()' function, which retrieves the the polygons that cross each other
# The parameter 'remove_self=TRUE' desconsider a polygon to be neighbor of itself.
intersect_list <- st_intersects(st_as_sf(nodes[1:9,]),remove_self=TRUE)
# st_intersects() returns a 'sgbp' list which contains the indices of the intersected 'sf' objects (neighbors)
intersect_list
# see the indeces of the neighbors of the first polygon:
intersect_list[[1]]
# [1] 6 (only one neighbor in this case, with index 6)

# We can convert the intersect list with a dataframe containing all the neighbor pairs
as.data.frame(intersect_list)
# 'row.id' contain the index of the reference 'sf' object and the 'col.id' have the index of its neighbor
################
# Lets visualize the neighbors obtained by 'st_intersects()'
# 'st_join()' returns a data frame with neighbor pairs from the two given 'sf' objects.
# In this case, the 'sf' objects are the same. The join condition is the 'st_intersects'.
neighbors_df <- st_join(st_as_sf(nodes[1:9,]), st_as_sf(nodes[1:9,]), join = st_intersects,remove_self=TRUE)
# The new 'status' column will only store if a given polygon have neighbor or not
neighbors_df$status <- ifelse(is.na(neighbors_df$vertex.y), "No", "Yes")

ggplot(data = neighbors_df)+
  geom_sf(aes(fill=status))+
  geom_text(aes(x=x.x,y=y.x,label = vertex.x))+
  scale_fill_manual(values = c("Yes"="#2b8cbe","No"="grey"),
                    na.value = "white")+
  theme_minimal()+
  labs(fill="Have Neighbor?",x="x",y="y")
# As you can see, not all cells were assigned to have neighbors because the their boundaries do not overlap. 
##########
# Expanding the polygons area to identify the neighboring cells with which the 'sf' polygon does not overlap
nodes <- st_as_sf(nodes[1:9,]) # just assigning the example dataframe as a 'sf' objects

# 'st_is_within_distance()' receives two 'sf' objects and retrieves the same 'sgbp' list returned by the 'st_intersects()' function.
# The 'dist' parameter is a distance used "expand" the polygons boundaries to assess if they would overlap.
# If yes, the indices of the overlaping polygons are stored.
neighbor_list <- st_is_within_distance(nodes,nodes, dist = 0.0005,remove_self=TRUE)
neighbor_list
# Here we can see that the neighbors were assigned properly, with the 'dist' equals to 0.0005.

# Inspecting the indices of the neighbors from the first cell/polygon.
neighbor_list[[1]]
# [1] 2 6 7
# Assessing the number of neighbors the first cell/polygons have.
length(neighbor_list[[1]])
# [1] 3

# Is possible to convert the convert the neighbors list to a dataframe:
neighbors_df <- as.data.frame(neighbor_list)
neighbors_df
########
# To illustrate the test made by 'st_is_within_distance()' we can use the 'st_buffer()' function.
# 'st_buffer()' change the polygon boundaries by a given distance ('dist') parameter
buffered_nodes <- st_buffer(st_as_sf(nodes[1:9,]),dist = 0.0005)
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
# Now all cells were assigned to have neighbors