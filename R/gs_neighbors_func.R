# 'gs_neighbors()' receives a GraphSpace object that do not has edges and has cell segmentation
# polygons in the nodes slot. The function also receives a dist parameter used defien a buffer region
# in the boundart of the cells.
gs_neighbors <- function(gs,dist){
  # Finding neighbors by cell segmentation
  # 'st_is_within_distance()' receives two 'sf' objects and retrieves the same 'sgbp' list returned by the 'st_intersects()' function.
  # The 'dist' parameter is a distance used "expand" the polygons boundaries to assess if they would overlap.
  # If yes, the indices of the overlaping polygons are stored.
  neighbor_list <- st_is_within_distance(st_as_sf(gs@nodes),
                                         st_as_sf(gs@nodes),
                                         dist = dist,
                                         remove_self=TRUE)
  # The 'neighbor_list' is a 'sgbp' list which contains the indices of the intersected 'sf' objects (neighbors)
  
  # Convert the intersect list with a dataframe containing all the neighbor pairs
  neighbor_df<- as.data.frame(neighbor_list)
  colnames(neighbor_df) <- c("from","to")
  
  neighbor_df <- neighbor_df %>%
    mutate(from_new = pmin(from, to),
      to_new   = pmax(from, to)) %>%
    select(-from, -to) %>%
    rename(from = from_new,to = to_new) %>%
    distinct()
  
  neighbor_df$from <- gs_vertex_attr(gs,"name")[neighbor_df$from]
  neighbor_df$to <- gs_vertex_attr(gs,"name")[neighbor_df$to]
  
  
  # Adding edges to the GraphSpace object by cell-cell proximity
  gs_neighbors <- gs |> gs_add_edges(neighbor_df)
  
  return(gs_neighbors)
}
