# Xenium Cell Com book

The complete source is **Xenium Cell Com book.qmd**. It produces one portable
HTML volume with five numbered chapters, in this order:

1. Inspecting neighbors — `InspectingNeighbors.Rmd`
2. Polygon geometry — `polygonGeometry.Rmd`
3. Delaunay and proximity networks — `DelaunayNet.Rmd`
4. Setting the proximity distance — `SettingDist_proximityNetwork.Rmd`
5. Visualizing direct neighbors — `Visualizing_DirectNeighbors.Rmd`

From this folder, render with:

```powershell
.\render-book.ps1
```

The helper finds Quarto (including RStudio's bundled copy) and uses a local
temporary folder to avoid Windows temporary-path issues. The source also works
with `quarto render "Xenium Cell Com book.qmd" --to html` in a configured R
and Quarto environment.

Use an R installation with `knitr` and `rmarkdown` available to Quarto. The
source checks its analysis dependencies and installs missing packages into
`Book/_R_library`. If Quarto is bundled with RStudio rather than on PATH, use
RStudio's Render button or the bundled Quarto executable.

The data file is read from `../data/Xenium_GS_Crop.RDS`; it is downloaded only
if missing. Generated figures, GIF frames, and analysis tables go into
`figures/`, `frames/`, and `results/`. The rendered HTML embeds its resources.
The ragg graphics device keeps figure generation fast on Windows.
System fonts avoid external font downloads. An executed Markdown intermediate
is retained to help diagnose any later formatting issues.

The book consolidates shared setup, computes spatial metrics and candidate
distances once, preserves raw-count correlations as exploratory analyses,
uses symmetric correlation color limits, applies one significance rule to
both correlation methods, compares edges by cell identity, and stores complete
plot objects in a list for the distance animation. Existing vignettes and
the project's original image files are not edited by the render.
