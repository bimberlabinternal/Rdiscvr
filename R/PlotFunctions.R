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
  return(FormatUmapForPresentation(P1) + scale_colour_gradientn(colours = c("navy", "dodgerblue", "gold", "red")))
}

#' @title theme_bimber
#'
#' @description Base theme for formatting ggplot objects using a standardized theme.
#'
#' @param standardPlot Logical. If TRUE, applies a cartesian article-style theme
#'   with axis titles/ticks/text; if FALSE, applies a minimal theme suited to
#'   donut/pie or other non-standard panels without axes.
#' @param legendPosition Character. Position of the legend; e.g., "right",
#'   "left", "top", "bottom", or "none".
#' @param face Character. Font face for titles, subtitles, and facet strips
#'   (e.g., "plain", "bold", "italic", "bold.italic").
#' @param titleJust Character. Title/subtitle horizontal justification:
#'   "left", "center", or "right".
#' @param baseSize Numeric. Base font size passed to the base theme.
#' @param axisTextsize Numeric. Font size for axis titles and tick labels when
#'   standardPlot = TRUE.
#' @param forceLegendAlpha Logical. If TRUE, forces legend keys to full opacity
#'   (override.aes alpha = 1).
#' @param adjustLegendKeySize Logical. If TRUE, enlarges legend key width/height
#'   for readability.
#'
#' @return ggplot object
#' @export

theme_bimber <- function(
    standardPlot = TRUE,
    legendPosition = "right",
    face = "bold",
    titleJust = "center",
    baseSize = 12,
    axisTextsize = 16,
    forceLegendAlpha = TRUE,
    adjustLegendKeySize = TRUE
){
  if(titleJust == "center"){
    hjust <-  .5
  } else if(titleJust == "right"){
    hjust <- 1
  } else if(titleJust == "left"){
    hjust <- 0
  }
  
  common_layer <- theme(
    text               = element_text(family = "Helvetica"),
    legend.position    = legendPosition,
    legend.title       = element_text(size = 16, face = face),
    legend.text        = element_text(size = 14),
    plot.title         = element_text(size = 22, face = face, hjust = hjust),
    plot.subtitle      = element_text(size = 18, face = face, hjust = hjust),
    strip.text.x       = element_text(size = 16, color = "black", face = face),
    strip.text.y       = element_text(size = 16, color = "black", face = face),
    strip.background   = element_rect(fill = "white", linewidth = 2)
  )
  
  if (standardPlot) {
    base_theme <- egg::theme_article(base_size = baseSize)
  } else {
    base_theme <- theme_void(base_size = baseSize)
  }
  
  if (!(standardPlot)) {
    axis_layer <- theme(
      axis.title  = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      axis.line   = element_blank(),
      panel.grid  = element_blank()
    )
  } else {
    axis_layer <- theme(
      axis.title  = element_text(size = axisTextsize),
      axis.text.x = element_text(size = axisTextsize, color = "black"),
      axis.text.y = element_text(size = axisTextsize, color = "black")
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
  
  list(base_theme, common_layer, axis_layer, guides_layer, legend_key_layer)
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
#' @param xLabel Character. Text label placed at the midpoint of the x-axis arrow.
#' @param yLabel Character. Text label placed at the midpoint of the y-axis arrow.
#' @param arrowLabelOffset Numeric. Offset of labels from the arrow line as a
#'   fraction of the respective axis range (vertical for x label, horizontal for y label).
#' @param arrowYlabelAngle Numeric. Rotation angle (degrees) for the y-axis arrow label.
#'
#' @return ggplot object
#' @export

ConvertAxesToArrows <- function(
    arrowLength = 0.3,
    arrowLinewidth = 1.5,
    arrowHeadLenMm = 2,
    arrowHeadAngle = 50,
    xLabel = "X_Axis",
    yLabel = "Y_Axis",
    arrowLabelOffset = 0.02,
    arrowYlabelAngle = 90
){
  structure(list(
    arrowLength = arrowLength,
    arrowLinewidth = arrowLinewidth,
    arrowHeadLenMm = arrowHeadLenMm,
    arrowHeadAngle = arrowHeadAngle,
    xLabel = xLabel,
    yLabel = yLabel,
    arrowLabelOffset = arrowLabelOffset,
    arrowYlabelAngle = arrowYlabelAngle
  ), class = "bimber_arrows")
}

#' @title ggplot_add method for bimber_arrows
#'
#' @description Internal helper used to add arrow-based axes and their labels to
#'   ggplot objects.
#'
#' @param object A code object created by {ConvertAxesToArrows()},containing arrow and label parameters.
#' @param plot A ggplot object to which the arrow-based axes and labels will be added.
#' @param object_name Character. The name of the {bimber_arrows} object
#'   being added (supplied by ggplot2; typically not used directly).
#'
#' @return A ggplot object with arrow-based axes and corresponding labels added.
#'
#' @method ggplot_add bimber_arrows
#' @export

ggplot_add.bimber_arrows <- function(object, plot, object_name){
  plot <- plot +
    coord_cartesian(clip = "off") +
    theme(axis.line = element_blank(),
          plot.margin = margin(6, 6, 16, 6))
  
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
  xLabel_y <- y0 - diff(yr) * object$arrowLabelOffset
  yLabel_x <- x0 - diff(xr) * object$arrowLabelOffset
  
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
      "text", x = xm, y = xLabel_y, label = object$xLabel,
      vjust = 1, hjust = 0.5
    ) +
    annotate(
      "text", x = yLabel_x, y = ym, label = object$yLabel,
      vjust = 0.5, hjust = 0.5, angle = object$arrowYlabelAngle
    )
}

#' @title ApplyBimberTheme
#' 
#' @description Wrapper for Rdiscvr plotting themes that formats ggplot objects using a standardized theme.
#' @param  plot a ggplot object
#' @param  standardPlot boolean controls theming defaults dependant on a plot is a UMAP, DonutPlot, or other non-Standard plot
#' @param legendPosition character, default `"right"`. Legend position passed to
#'   `theme(legend.position = ...)`. Supported values: `"none"`, `"left"`, `"right"`,
#'   `"bottom"`, `"top"`.
#' @param face character, default `"bold"`. Font face for plot titles, subtitles,
#'   and facet strips (e.g., `"plain"`, `"bold"`, `"italic"`, `"bold.italic"`).
#' @param titleJust character, default `"center"`. Horizontal justification of
#'   title/subtitle. One of `"left"`, `"center"`, `"right"`.
#' @param baseSize numeric, default `12`. Base font size fed to the base theme.
#' @param axisTextsize numeric, default `16`. Font size for axis titles and tick
#'   labels when `standardPlot = TRUE`. Also used to scale arrow-label text.
#' @param forceLegendAlpha logical, default `TRUE`. If `TRUE`, legend guide
#'   overrides mapped alpha to `1` so legend keys are fully opaque.
#' @param adjustLegendKeySize logical, default `TRUE`. If `TRUE`, increases
#'   legend key width/height for readability.
#' @param useArrows logical, default `FALSE`. If `TRUE`, removes default axis
#'   lines and draws perpendicular x/y arrows at the lower-left
#'   of the panel with customizable head/shaft and centered labels.
#'
#' @param arrowLength numeric in (0,1), default `0.3`. Arrow length as
#'   a fraction of the respective axis range (x and y independently).
#' @param arrowLinewidth numeric, default `1.5`. Line width of the arrow shafts.
#' @param arrowHeadLenMm numeric, default `2`. Arrowhead length in millimeters.
#' @param arrowHeadAngle numeric, default `50`. Arrowhead angle (larger =
#'   “bigger” head).
#'
#' @param xLabel character, default `"X_Axis"`. Text shown at the midpoint of
#'   the x-axis arrow.
#' @param yLabel character, default `"Y_Axis"`. Text shown at the midpoint of
#'   the y-axis arrow.
#' @param arrowLabelOffset numeric, default `0.02`. Offset of the label text
#'   from the arrow shafts, expressed as a fraction of the corresponding axis
#'   range (x-label offset is vertical; y-label offset is horizontal).
#' @param arrowYlabelAngle numeric, default `90`. Rotation for the y-axis arrow
#'   label (set to `90` for vertical).
#' @return ggplot object
#' @export

ApplyBimberTheme <- function(plot = NULL,
                             standardPlot = TRUE,
                             legendPosition = "right",
                             face = "bold",
                             titleJust = "center",
                             baseSize = 12,
                             axisTextsize = 16,
                             forceLegendAlpha = TRUE,
                             adjustLegendKeySize = TRUE,
                             useArrows = FALSE,
                             arrowLength = 0.3,
                             arrowLinewidth = 1.5,
                             arrowHeadLenMm = 2,
                             arrowHeadAngle = 50,
                             xLabel = "X_Axis",
                             yLabel = "Y_Axis",
                             arrowLabelOffset = 0.02,
                             arrowYlabelAngle = 90){
  
  #################
  ### Sanitize  ###
  #################
  
  if(is.null(plot)){
    stop("No ggplot object supplied, please provide a ggplot object.")
  } else if (!is.logical(standardPlot)){
    stop("standardPlot argument is not logical. Please set to TRUE/FALSE")
  }
  
  plot <- plot +
    theme_bimber(standardPlot = standardPlot,
                 legendPosition = legendPosition,
                 face = face,
                 titleJust = titleJust,
                 baseSize = baseSize,
                 axisTextsize = axisTextsize,
                 forceLegendAlpha = forceLegendAlpha,
                 adjustLegendKeySize = adjustLegendKeySize)
  
  if (useArrows) {
    plot <- plot + ConvertAxesToArrows(
      arrowLength = arrowLength,
      arrowLinewidth = arrowLinewidth,
      arrowHeadLenMm = arrowHeadLenMm,
      arrowHeadAngle = arrowHeadAngle,
      xLabel = xLabel,
      yLabel = yLabel,
      arrowLabelOffset = arrowLabelOffset,
      arrowYlabelAngle = arrowYlabelAngle
    )
  }
  
  return(plot)
}
