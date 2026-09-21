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
dist_vec <- seq(from = 0, to = 0.002, by = 0.0001)
edges_by_dist <- data.frame(
  dist = dist_vec,
  edges = NA,
  row.names = dist_vec
)

g <- ggplot(gs_crop) +
  geom_edgespace(color="red")+
  geom_nodespace(color="black",
                 size = .2,pch = 19) +
  scale_colour_continuous(palette = cpal, limits = data_range)+
  theme_gspace_coords(theme = "th3",is_norm = TRUE)+
  theme_minimal()+
  labs(title = "Cell proximity, No edges",x=NULL,y=NULL)
png_files <- c(g)

inicio <- Sys.time()
for (i in 1:length(dist_vec)) {
  gs_prox <- gs_neighbors(gs_crop,dist=dist_vec[i])
  
  g <- ggplot(gs_prox) +
    geom_edgespace(color="red")+
    geom_nodespace(color="black",
                   size = .2,pch = 19) +
    scale_colour_continuous(palette = cpal, limits = data_range)+
    theme_gspace_coords(theme = "th3",is_norm = TRUE)+
    theme_minimal()+
    labs(title = paste("Cell proximity, dist=",dist_vec[i]),x=NULL,y=NULL)
  png_files <- c(png_files,g)
  
  edges_by_dist[i,]$edges <-ecount(gs_prox@graph)
  print(paste("Iteration:",i,"Dist = ",dist_vec[i],"OK"))
}
fim <- Sys.time()
tempo <- fim - inicio
tempo

file_paths <- sapply(seq_along(png_files), function(i) {
  file_path <- file.path("images/frames/", sprintf("frame_%03d.png", i))
  
  # ggsave handles saving the ggplot object automatically
  ggsave(
    filename = file_path, 
    plot = png_files[[i]], 
    width = 6,      # Width in inches
    height = 6,     # Height in inches
    dpi = 100       # 6 inches * 100 dpi = 600x600 pixels
  )
  
  return(file_path)
})

gifski(png_files = file_paths[-1],
       gif_file = "images/CellProx_Dist.gif",
       width = 600,
       height = 600,
       delay = 0.2)

g<- ggplot(edges_by_dist,aes(x=dist,y=edges))+
         geom_line()+
         geom_point()+
  labs(y="Number of edges")

ggsave("images/Proximity_Dist_numberEdges.png",g,
       height = 4,width = 7,units = "in",dpi = 120)
######################################################
# Calculate the greates distance btween any pair of its vertices
# load function that calculate the maximum distance of a polygon (polygon diameter)
source("R/maxDistPol_func.R")

# Assign maximum distance (polygon diameter)
all_geom <- gs_crop@nodes %>% select(vertex,geometry)

for(i in 1:length(all_geom$geometry)){
  all_geom$maxDist[i] <- maxDistPol(all_geom$geometry[i])
}
gs_vertex_attr(gs_crop, "maxDist") <- all_geom$maxDist
rm(all_geom)

# Store the meanMaxDist between all the cell polygons
meanMaxDist <- mean(gs_vertex_attr(gs_crop, "maxDist"))
######################################################        
0.0005/meanMaxDist
