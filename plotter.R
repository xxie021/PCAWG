library(ggplot2)

# Constants
kPlotOutPath.default <- "plot_out/"

kPlotTheme.default <- theme(plot.title = element_text(face = "bold", size = 16,
                                                      hjust = 0.5,
                                                      margin = margin(b = 15)),
                            plot.margin = unit(rep(15, 4), units = "pt"))

kBarTheme <- theme(
  axis.title.y = element_blank(),
  axis.text = element_text(face = "bold", colour = "black", size = 11)
)

kBarTheme.num <- theme(
  axis.title.x = element_text(face = "bold", size = 12, margin = margin(15))
)

kBarTheme.pct <- theme(
  axis.title.x = element_blank()
)

kStackedBarTheme <- theme(
  axis.title.x = element_text(face = "bold", size = 14, margin = margin(15)),
  axis.title.y = element_blank(),
  axis.text = element_text(face = "bold", colour = "black"),
  axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
  axis.text.y = element_text(size = 11),
  legend.background = element_rect(colour = "gray"),
  legend.title = element_text(face = "bold", size = 10.5),
  legend.title.align = 0.5
)

kHeatmapTheme <- theme(
  axis.title = element_text(face = "bold", colour = "black", size = 12),
  axis.title.x = element_text(margin = margin(t = 12)),
  axis.title.y = element_text(margin = margin(r = 12)),
  axis.text = element_text(face = "bold", colour = "black", size = 8),
  legend.background = element_rect(colour = "gray"),
  legend.margin = margin(5, 25, 5, 25),
  legend.title = element_text(face = "bold", size = 10.5),
  legend.title.align = 0.5,
  legend.position = "bottom",
  panel.spacing = unit(0, "pt"),
  strip.background = element_rect(fill = "black", colour = "goldenrod",
                                  size = 1),
  strip.text = element_text(face = "bold", colour = "white", size = 9)
)

# Inheritance
# Top: XPlotter
# 2nd Level: PdfPlotter, JpgPlotter
XPlotter <- function(plotter.type, plot.out.path = kPlotOutPath.default,
                     plot.theme.extra = NA) {
  me <- list(
    type = plotter.type,
    out.path = plot.out.path
  )
  
  if (is.theme(plot.theme.extra)) {
    me$theme <- kPlotTheme.default + plot.theme.extra
  } else {
    me$theme <- kPlotTheme.default
  }
  
  class(me) <- append(class(me), "XPlotter")
  return(me)
}

PdfPlotter <- function(plot.out.path = kPlotOutPath.default,
                       plot.theme.extra = NA) {
  me <- XPlotter("pdf", plot.out.path, plot.theme.extra)
  
  class(me) <- append(class(me), "PdfPlotter")
  return(me)
}

JpgPlotter <- function(plot.out.path = kPlotOutPath.default,
                       plot.theme.extra = NA) {
  me <- XPlotter("jpg", plot.out.path, plot.theme.extra)
  
  class(me) <- append(class(me), "JpgPlotter")
  return(me)
}

Summary <- function(obj, ..., out.path = NA) {
  if (!is.na(out.path) && !dir.exists(out.path)) {
    dir.create(out.path)
  }
  UseMethod("Summary", obj)
}

PlotBaseSpectrum <- function(obj, plotter, ...) {
  if (class(plotter)[2] != "XPlotter") {
    stop("Invaild XPlotter found", call. = F)
  }
  
  if (!dir.exists(plotter$out.path)) {
    dir.create(plotter$out.path)
  }
  UseMethod("PlotBaseSpectrum", obj)
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
