#' @title FormatUmapForPresentation
#'
#' @description Calls renames the axes to UMAP_1 and UMAP_2 on a DimPlot/FeaturePlot
#' @param P1 The plot object
#' @return A plot object
#' @export
FormatUmapForPresentation <- function(P1) {
  return(P1 & labs(x = 'UMAP_1', y = 'UMAP_2'))
}

#' @title FormatFeaturePlotColorScale
#'
#' @description Calls FormatFeaturePlotColorScale and sets a navy/red color scale for Seurat FeaturePlot objects
#' @param P1 The plot object
#' @return A plot object
#' @export
FormatFeaturePlotColorScale <- function(P1) {
  return(FormatUmapForPresentation(P1) + ggplot2::scale_colour_gradientn(colours = c("navy", "dodgerblue", "gold", "red")))
}