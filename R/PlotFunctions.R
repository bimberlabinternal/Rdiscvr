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

#' @title ApplyBimberTheme
#' 
#' @description Formats ggplot objects using a standardized theme.
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
#' @param arrowLength numeric in (0,1), default `0.3`. Arrow **shaft length** as
#'   a fraction of the respective axis range (x and y independently).
#' @param arrowLinewidth numeric, default `1.5`. Line width of the arrow shafts.
#' @param arrowHeadLenMm numeric, default `2`. Arrowhead length in millimeters.
#' @param arrowHeadAngle numeric, default `50`. Arrowhead angle (larger =
#'   “fatter” head).
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
#' 
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
                             # --- arrow controls ---
                             arrowLength = 0.3,          
                             arrowLinewidth = 1.5,      
                             arrowHeadLenMm = 2,      
                             arrowHeadAngle = 50,      
                             # --- NEW: labels on arrows ---
                             xLabel = "X_Axis",
                             yLabel = "Y_Axis",
                             arrowLabelOffset = 0.02,  
                             arrowYlabelAngle = 90) {
  
  #################
  ### Sanitize  ###
  #################
  
  if(is.null(plot)){
    stop("No ggplot object supplied, please provide a ggplot object.")
  } else if (!is.logical(standardPlot)){
    stop("standardPlot argument is not logical. Please set to TRUE/FALSE")
  }
  
  
  if(titleJust == "center"){
    hjust <-  .5
  } else if(titleJust == "right"){
    hjust <- 1
  } else if(titleJust == "left"){
    hjust <- 0
  }
  
  #################
  ### BasePlot Layer  ###
  #################
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
  
  #################
  ### Options  ###
  #################
  if (standardPlot) {
    base_theme <- egg::theme_article(baseSize = baseSize)
  } else {
    base_theme <- theme_void(baseSize = baseSize)
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
    guides_layer <- ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  } else {
    guides_layer <- NULL
  }
  
   if (adjustLegendKeySize) {
    legend_key_layer <- ggplot2::theme(
      legend.key.height = grid::unit(10, "pt"),
      legend.key.width  = grid::unit(10, "pt")
    )
   } else {
     legend_key_layer <- NULL
   }
  
  plot <- plot + base_theme + common_layer + axis_layer + guides_layer + legend_key_layer
  
  if (useArrows) {
    # allow drawing outside panel; hide default axis lines
    plot <- plot +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(axis.line = ggplot2::element_blank())
    
    # build to get panel ranges
    gb <- ggplot2::ggplot_build(plot)
    pp <- gb$layout$panel_params[[1]]
    
    xr <- if (!is.null(pp$x.range)) pp$x.range else pp$x$range$range
    yr <- if (!is.null(pp$y.range)) pp$y.range else pp$y$range$range
    
    x0 <- xr[1]; y0 <- yr[1]
    dx <- diff(xr) * arrowLength
    dy <- diff(yr) * (arrowLength + 0.05)
    
    # midpoints of the shafts
    
    xm <- x0 + dx/2
    ym <- y0 + dy/2
    
    # put X label *below* the x-axis arrow (negative offset from y0)
    xLabel_y <- y0 - diff(yr) * arrowLabelOffset   
    yLabel_x <- x0 - diff(xr) * arrowLabelOffset
    
    # make sure there is room below for the text
    plot <- plot + ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 16, 6))
    
    arrow_spec <- grid::arrow(type = "closed", ends = "last",
                              length = grid::unit(arrowHeadLenMm, "mm"),
                              angle = arrowHeadAngle)
    
    plot <- plot +
      # shafts with heads
      ggplot2::annotate("segment",
                        x = x0, y = y0, xend = x0 + dx, yend = y0,
                        linewidth = arrowLinewidth, arrow = arrow_spec
      ) +
      ggplot2::annotate("segment",
                        x = x0, y = y0, xend = x0, yend = y0 + dy,
                        linewidth = arrowLinewidth, arrow = arrow_spec
      ) +
      # centered labels
      ggplot2::annotate(
        "text", x = xm, y = xLabel_y, label = xLabel,
        size = axisTextsize/3, vjust = 1, hjust = 0.5
      ) +
      ggplot2::annotate(
        "text", x = yLabel_x, y = ym, label = yLabel,
        size = axisTextsize/3, vjust = 0.5, hjust = 0.5, angle = arrowYlabelAngle
      )
  }
  
  return(plot)
}