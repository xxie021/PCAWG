library(ggplot2)

# Constants
kPlotOutPath.default <- "plot_out/"

kPlotTheme.default <- theme(plot.title = element_text(face = "bold", size = 16,
                                                      hjust = 0.5,
                                                      margin = margin(b = 15)),
                            plot.margin = unit(rep(15, 4), units = "pt"))
kAxisTheme.default <- theme(axis.text = element_text(face = "bold",
                                                     colour = "black"),
                            axis.text.y = element_text(size = 11),
                            axis.title.y = element_blank())
kLegendTheme.default <- theme(legend.background = element_rect(colour = "gray"),
                              legend.title = element_text(face = "bold",
                                                          size = 10.5),
                              legend.title.align = 0.5)

# Inheritance
# Top: XPlotter
# 2nd Level: PdfPlotter, JpgPlotter
XPlotter <- function(plotter.type, plotter.out.path = kPlotOutPath.default,
                     plot.theme = kPlotTheme.default, plot.theme.extra = NA,
                     axis.theme = kAxisTheme.default, axis.theme.extra = NA,
                     legend.theme = kLegendTheme.default,
                     legend.theme.extra = NA) {
  if (is.theme(plot.theme) && is.theme(plot.theme.extra)) {
    plot.theme <- plot.theme + plot.theme.extra
  } else if (!is.theme(plot.theme) && is.theme(plot.theme.extra)) {
    plot.theme <- plot.theme.extra
  } else if (!is.theme(plot.theme) && !is.theme(plot.theme.extra)) {
    plot.theme <- NA
  }
  
  if (is.theme(axis.theme) && is.theme(axis.theme.extra)) {
    axis.theme <- axis.theme + axis.theme.extra
  } else if (!is.theme(axis.theme) && is.theme(axis.theme.extra)) {
    axis.theme <- axis.theme.extra
  } else if (!is.theme(axis.theme) && !is.theme(axis.theme.extra)) {
    axis.theme <- NA
  }
  
  if (is.theme(legend.theme) && is.theme(legend.theme.extra)) {
    legend.theme <- legend.theme + legend.theme.extra
  } else if (!is.theme(legend.theme) && is.theme(legend.theme.extra)) {
    legend.theme <- legend.theme.extra
  } else if (!is.theme(legend.theme) && !is.theme(legend.theme.extra)) {
    legend.theme <- NA
  }
  
  me <- list(
    type = plotter.type,
    out.path = plotter.out.path,
    plot.theme = plot.theme,
    axis.theme = axis.theme,
    legend.theme = legend.theme
  )
  
  class(me) <- append(class(me), "XPlotter")
  return(me)
}

PdfPlotter <- function(plotter.out.path = kPlotOutPath.default,
                       plot.theme = kPlotTheme.default, plot.theme.extra = NA,
                       axis.theme = kAxisTheme.default, axis.theme.extra = NA,
                       legend.theme = kLegendTheme.default,
                       legend.theme.extra = NA) {
  me <- XPlotter("pdf", plotter.out.path, plot.theme, plot.theme.extra,
                 axis.theme, axis.theme.extra, legend.theme, legend.theme.extra)
  
  class(me) <- append(class(me), "PdfPlotter")
  return(me)
}

JpgPlotter <- function(plotter.out.path = kPlotOutPath.default,
                       plot.theme = kPlotTheme.default,
                       plot.theme.extra = NA,
                       axis.theme = kAxisTheme.default,
                       axis.theme.extra = NA,
                       legend.theme = kLegendTheme.default,
                       legend.theme.extra = NA) {
  me <- XPlotter("jpg", plotter.out.path, plot.theme, plot.theme.extra,
                 axis.theme, axis.theme.extra, legend.theme, legend.theme.extra)
  
  class(me) <- append(class(me), "JpgPlotter")
  return(me)
}

PlotBarData <- function(obj, plotter, ...) {
  if (class(plotter)[2] != "XPlotter") {
    stop("Invaild XPlotter found", call. = F)
  }
  
  if (!dir.exists(plotter$out.path)) {
    dir.create(plotter$out.path)
  }
  UseMethod("PlotBarData", obj)
}

PlotHeatmap <- function(obj, plotter, ...) {
  if (class(plotter)[2] != "XPlotter") {
    stop("Invaild XPlotter found", call. = F)
  }
  
  if (!dir.exists(plotter$out.path)) {
    dir.create(plotter$out.path)
  }
  UseMethod("PlotHeatmap", obj)
}
