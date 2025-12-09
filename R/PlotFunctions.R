#' @title RenameUmapAxes
#'
#' @param prefix The prefix to apply to the axis labels
#' @description Renames the axes to UMAP_1 and UMAP_2 on a DimPlot/FeaturePlot. Usage: Dimplot() + RenameUmapAxes()
#' @return A plot object
#' @export
#'
#' @examples
#' \dontrun{
#' Seurat::DimPlot(seuratObj) +
#'      RenameUmapAxes()
#'
#' Seurat::DimPlot(seuratObj) +
#'      RenameUmapAxes(prefix = 'tSNE')
#' }

RenameUmapAxes <- function(prefix = 'UMAP') {
  return(ggplot2::labs(x = paste0(prefix, '_1'), y = paste0(prefix, '_2')))
}

#' @title FormatFeaturePlotColorScale
#'
#' @description Calls FormatFeaturePlotColorScale and sets a navy/red color scale for Seurat FeaturePlot objects
#' @param prefix The string prefix for the axes. Passed to RenameUmapAxes()
#' @return A plot object
#' @export
#' @examples
#' \dontrun{
#' Seurat::DimPlot(seuratObj) +
#'      FormatFeaturePlotColorScale()
#' }
FormatFeaturePlotColorScale <- function(prefix = 'UMAP') {
  return(list(
           RenameUmapAxes(prefix) &
           scale_colour_gradientn(colours = c("navy", "dodgerblue", "gold", "red"))
  ))
}

#' @title theme_bimberlab
#'
#' @description Base theme for formatting ggplot objects using a standardized theme.
#'
#' @param minimalPlot Logical. If TRUE, theme_minimal() is applied and many plot elements are hidden.
#'   with axis titles/ticks/text; if FALSE, applies a minimal theme suited to
#'   donut/pie or other non-standard panels without axes.
#' @param legendPosition Character. Position of the legend; e.g., "right",
#'   "left", "top", "bottom", or "none".
#' @param fontFamily The default font family for this plot
#' @param face Character. Font face for titles, subtitles, and facet strips
#'   (e.g., "plain", "bold", "italic", "bold.italic").
#' @param titleJust Character. Title/subtitle horizontal justification:
#'   "left", "center", or "right".
#' @param baseSize Numeric. Base font size passed to the base theme.
#' @param axisTextSize Numeric. Font size for axis titles and tick labels.
#' @param forceLegendAlpha Logical. If TRUE, forces legend keys to full opacity
#'   (override.aes alpha = 1).
#' @param adjustLegendKeySize Logical. If TRUE, enlarges legend key width/height
#'   for readability.
#' @param ... Additional parameters that will be passed to theme()
#'
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' Seurat::DimPlot(seuratObj) +
#'      theme_bimberlab()
#'
#' # Pass additional theme() arguments like this:
#' Seurat::DimPlot(seuratObj) +
#'      theme_bimberlab(
#'          legend.position = 'none'
#'      )
#' }
theme_bimberlab <- function(
    minimalPlot = FALSE,
    legendPosition = "right",
    fontFamily = "Helvetica",
    face = "bold",
    titleJust = "center",
    baseSize = 12,
    axisTextSize = 16,
    forceLegendAlpha = TRUE,
    adjustLegendKeySize = TRUE,
    ...
){

  hjust <- switch(titleJust,
                  "center" = 0.5,
                  "right" = 1,
                  "left" = 0,
                  stop(paste0('Unknown value for titleJust: ', titleJust))
  )
  
  common_layer <- theme(
    # TODO: on windows this gives the error: Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : font family not found in Windows font database
    #text               = element_text(family = fontFamily),
    legend.position    = legendPosition,
    legend.title       = element_text(face = face),
    legend.text        = element_text(),
    plot.title         = element_text(face = face, hjust = hjust),
    plot.subtitle      = element_text(face = face, hjust = hjust),
    strip.text.x       = element_text(color = "black", face = face),
    strip.text.y       = element_text(color = "black", face = face),
    strip.background   = element_rect(linewidth = 2)
  )
  
  if (minimalPlot) {
    base_theme <- theme_void(base_size = baseSize)
    axis_layer <- theme(
      axis.title  = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      axis.line   = element_blank(),
      panel.grid  = element_blank()
    )
  } else {
    base_theme <- egg::theme_article(base_size = baseSize)
    axis_layer <- theme(
      axis.title  = element_text(size = axisTextSize),
      axis.text.x = element_text(size = axisTextSize, color = "black"),
      axis.text.y = element_text(size = axisTextSize, color = "black")
    )
  }
  
  if (forceLegendAlpha) {
    guides_layer <- guides(color = guide_legend(override.aes = list(size = 6, alpha = 1)))
  } else {
    guides_layer <- NULL
  }
  
  if (adjustLegendKeySize) {
    legend_key_layer <- theme(
      legend.key.height = grid::unit(10, "pt"),
      legend.key.width  = grid::unit(10, "pt")
    )
  } else {
    legend_key_layer <- NULL
  }
  
  list(base_theme, common_layer, axis_layer, guides_layer, legend_key_layer, theme(...))
}

#' @title ConvertAxesToArrows
#'
#' @description Replaces standard x/y axis lines with short, arrow-based axes and
#'   adds labels at their midpoints (useful for UMAPs and other DR plots).
#'
#' @param arrowLength Numeric. Shaft length as a fraction of the axis range.
#' @param arrowLinewidth Numeric. Line width of the arrow shafts.
#' @param arrowHeadLenMm Numeric. Arrowhead length in millimeters.
#' @param arrowHeadAngle Numeric. Arrowhead angle in degrees (larger = “bigger” head).
#' @param arrowLabelOffsetX Numeric. Offset of labels from the arrow line as a fraction of the x axis range (vertical for x label, horizontal for y label).
#' @param arrowLabelOffsetY Numeric. Offset of labels from the arrow line as a fraction of the x axis range (vertical for x label, horizontal for y label).
#' @param arrowYlabelAngle Numeric. Rotation angle (degrees) for the y-axis arrow label.
#' @param renameUmapAxes If true, RenameUmapAxes() will be run
#' @param plotAxisPrefix Passed directly to RenameUmapAxes(prefix = plotAxisPrefix)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' Seurat::DimPlot(seuratObj) +
#'      ConvertAxesToArrows()
#'
#' Seurat::DimPlot(seuratObj) +
#'      ConvertAxesToArrows(plotAxisPrefix = 'tSNE')
#' }
ConvertAxesToArrows <- function(
    arrowLength = 0.3,
    arrowLinewidth = 1.5,
    arrowHeadLenMm = 2,
    arrowHeadAngle = 50,
    arrowLabelOffsetX = 0.08,
    arrowLabelOffsetY = 0.02,
    arrowYlabelAngle = 90,
    renameUmapAxes = TRUE,
    plotAxisPrefix = 'UMAP'
){
  structure(list(
    arrowLength = arrowLength,
    arrowLinewidth = arrowLinewidth,
    arrowHeadLenMm = arrowHeadLenMm,
    arrowHeadAngle = arrowHeadAngle,
    arrowLabelOffsetX = arrowLabelOffsetX,
    arrowLabelOffsetY = arrowLabelOffsetY,
    arrowYlabelAngle = arrowYlabelAngle,
    renameUmapAxes = renameUmapAxes,
    plotAxisPrefix = plotAxisPrefix
  ), class = "dimplot_arrows")
}

#' @title ggplot_add method for dimplot_arrows
#'
#' @description Internal helper used to add arrow-based axes and their labels to
#'   ggplot objects.
#'
#' @param object A code object created by \code{ConvertAxesToArrows()},containing arrow and label parameters.
#' @param plot A ggplot object to which the arrow-based axes and labels will be added.
#' @param object_name Character. The name of the \code{dimplot_arrows} object
#'   being added (supplied by ggplot2; typically not used directly).
#'
#' @return A ggplot object with arrow-based axes and corresponding labels added.
#'
#' @method ggplot_add dimplot_arrows
#' @export

ggplot_add.dimplot_arrows <- function(object, plot, object_name){
  plot <- plot +
    coord_cartesian(clip = "off") +
    theme(axis.line = element_blank(),
          plot.margin = margin(6, 6, 16, 6),
          axis.title  = element_blank(),
          axis.text   = element_blank(),
          axis.ticks  = element_blank(),
          panel.grid  = element_blank()
    )

  if (object$renameUmapAxes) {
    plot <- plot + RenameUmapAxes(prefix = object$plotAxisPrefix)
  }

  gb <- ggplot_build(plot)
  pp <- gb$layout$panel_params[[1]]
  
  xr <- if (!is.null(pp$x.range)) pp$x.range else pp$x$range$range
  yr <- if (!is.null(pp$y.range)) pp$y.range else pp$y$range$range
  
  x0 <- xr[1]; y0 <- yr[1]
  dx <- diff(xr) * object$arrowLength
  dy <- diff(yr) * (object$arrowLength + 0.05)
  
  xm <- x0 + dx/2
  ym <- y0 + dy/2
  
  # below x-axis arrow, left of y-axis arrow
  xLabel_y <- y0 - diff(yr) * object$arrowLabelOffsetY
  yLabel_x <- x0 - diff(xr) * object$arrowLabelOffsetX
  
  arrow_spec <- grid::arrow(
    type = "closed", ends = "last",
    length = grid::unit(object$arrowHeadLenMm, "mm"),
    angle = object$arrowHeadAngle
  )
  
  plot +
    annotate("segment",
             x = x0, y = y0, xend = x0 + dx, yend = y0,
             linewidth = object$arrowLinewidth, arrow = arrow_spec
    ) +
    annotate("segment",
             x = x0, y = y0, xend = x0, yend = y0 + dy,
             linewidth = object$arrowLinewidth, arrow = arrow_spec
    ) +
    annotate(
      "text", x = xm, y = xLabel_y, label = ggplot2::get_labs(plot)$x,
      vjust = 1, hjust = 0.5
    ) +
    annotate(
      "text", x = yLabel_x, y = ym, label = ggplot2::get_labs(plot)$y,
      vjust = 0.5, hjust = 0.5, angle = object$arrowYlabelAngle
    )
}