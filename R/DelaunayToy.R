# 1. Install and load the required libraries
#install.packages(c("deldir", "igraph"))
library(deldir)
library(igraph)

# 2. Create sample 2D points (X and Y coordinates)
set.seed(42)
points <- data.frame(
  x = runif(20, 0, 100),
  y = runif(20, 0, 100)
)

# 3. Compute the Delaunay Triangulation
triangulation <- deldir(points$x, points$y)

# 4. Extract the edges (connections between points)
edges_matrix <- triangulation$delsgs[, c("ind1", "ind2")]

# 5. Build the igraph object
delaunay_graph <- graph_from_edgelist(as.matrix(edges_matrix), directed = FALSE)

# 6. Plot the resulting graph
plot(delaunay_graph, layout = as.matrix(points), 
     vertex.size = 7, vertex.label = NA, vertex.color = "tomato")
