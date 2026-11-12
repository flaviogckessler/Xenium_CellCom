library("RGraphSpace")
library("sf")
library("ggplot2")

gs_crop <- readRDS("data/gs_neig2.rds")
head(gs_crop[['edges']])

# Plot three built-in prototypes side by side
plot(GlyphArrow, GlyphCross1, GlyphCross2, GlyphChevron1, GlyphChevron2)

class(GlyphArrow)

ggplot(gs_crop) +
  geom_sf(mapping = aes(geometry = geometry,fill = log2(Slc17a7 + 1))) +
  #scale_fill_continuous(palette = cpal, limits = data_range) +
  geom_sf(data=gs_crop[["edges"]]$intersect_point,color="black",alpha=1,shape=8)+
  #geom_nodespace(GlyphArrow)+
  #my_theme+
  labs(title = "Point of physical contact between neighboring cells")
