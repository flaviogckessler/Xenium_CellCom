# Function: calculate the max distance polygon diameter
maxDistPol <- function(polygon){
  # Extract polygon points
  polygon_points <- st_cast(polygon, "POINT")
  polygon_points <- polygon_points[-1]
  
  # Find the absolute maximum distance between any two vertices
  max_distance <- max(st_distance(polygon_points, polygon_points))
  
  return(max_distance)
}