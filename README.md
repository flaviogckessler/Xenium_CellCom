# Xenium_CellCom

R workflows for exploring cell relationships in Xenium spatial transcriptomics
data, supporting cell–cell communication analysis. The examples use RGraphSpace
and cell segmentation polygons to inspect neighbors, measure cell geometry,
and compare Delaunay and proximity networks.

## Analysis book

[Read Xenium Cell-Cell Communication Book online](https://flaviogckessler.github.io/Xenium_CellCom/Book/XeniumCellCombook.html).
Its five chapters cover neighbor inspection, polygon geometry,
Delaunay networks, proximity distance thresholds, and direct neighbor
visualization. Figures and animation are embedded for offline reading.

See the [Quarto source](Book/XeniumCellCombook.qmd) and
[rendering instructions](Book/README.md) to reproduce the book.

## Repository contents

- [R](R/): analysis scripts and helper functions.
- [Vignettes](Vignettes/): R Markdown walkthroughs of the scripts.
- [Book](Book/): the integrated Quarto source, HTML book, and render helper.
- [data](data/) and [images](images/): example data and analysis figures.
