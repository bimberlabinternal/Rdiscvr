library(testthat)
library(Seurat)

test_that("Plot extensions do not error", {
  seuratObj <- readRDS('../testdata/small_RIRA.rds')

  ThePlot <- Seurat::DimPlot(seuratObj)

  ThePlot +
    theme_bimberlab()

  ThePlot +
    RenameUmapAxes()
  
  ThePlot +
    RenameUmapAxes() +
    theme_bimberlab()

  ThePlot +
    theme_bimberlab() +
    RenameUmapAxes()
  
  ThePlot +
    ConvertAxesToArrows()

  ThePlot +
    ConvertAxesToArrows(arrowLabelOffsetX = 0.08, arrowLabelOffsetY = 0.02) +
    theme_bimberlab()

  FP <- FeaturePlot(seuratObj, features = 'IFNG')

  FP +
    theme_bimberlab(
      legend.position = 'none'
    )

  FP +
    ConvertAxesToArrows()

  FP +
    FormatFeaturePlotColorScale()

  FP +
    ConvertAxesToArrows() +
    theme_bimberlab() +
    FormatFeaturePlotColorScale()

  FP +
    FormatFeaturePlotColorScale() +
    ConvertAxesToArrows() +
    theme_bimberlab()
})

test_that("Plot extensions work for simple ggplots", {
  P <- ggplot(datasets::mtcars, aes(x = mpg, y = cyl)) +
    geom_point()
  
  P
  P +
    ConvertAxesToArrows()
  
  P +
    theme_bimberlab()
  
  P +
    theme_bimberlab() +
    RenameUmapAxes()
  
})