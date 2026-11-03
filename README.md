# Xenium_CellCom

R workflows for exploring cell relationships in Xenium spatial transcriptomics
data, supporting cell–cell communication analysis. The examples use RGraphSpace
and cell segmentation polygons to inspect neighbors, measure cell geometry,
and compare Delaunay and proximity networks.

## Analysis book

[Download the HTML book](Book/Xenium%20Cell%20Com%20book.html) and open it in a
browser. Its five chapters cover neighbor inspection, polygon geometry,
Delaunay networks, proximity distance thresholds, and direct neighbor
visualization. Figures and animation are embedded for offline reading.

See the [Quarto source](Book/Xenium%20Cell%20Com%20book.qmd) and
[rendering instructions](Book/README.md) to reproduce the book.

## Repository contents

- [R](R/): analysis scripts and helper functions.
- [Vignettes](Vignettes/): R Markdown walkthroughs of the scripts.
- [Book](Book/): the integrated Quarto source, HTML book, and render helper.
- [data](data/) and [images](images/): example data and analysis figures.
