library(ggplot2)

# Constants
kSupportedOutputFormats <- c("PdfFileOut", "JpgFileOut")

kDefaultTheme <- theme(plot.title = element_text(face = "bold", size = 16,
                                                 hjust = 0.5,
                                                 margin = margin(b = 15)),
                       plot.margin = unit(rep(15, 4), units = "pt"))

kBarTheme.default <- theme(
  axis.title.y = element_blank(),
  axis.text = element_text(face = "bold", colour = "black", size = 11)
)

kBarTheme.num <- theme(
  axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 15))
)

kBarTheme.pct <- theme(
  axis.title.x = element_blank()
)

kStackedBarTheme.default <- theme(
  axis.title.x = element_text(face = "bold", size = 14,
                              margin = margin(t = 15)),
  axis.title.y = element_blank(),
  axis.text = element_text(face = "bold", colour = "black"),
  axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
  axis.text.y = element_text(size = 11)
)

kHeatmapTheme.default <- theme(
  axis.title = element_text(face = "bold", colour = "black", size = 12),
  axis.title.x = element_text(margin = margin(t = 12)),
  axis.title.y = element_text(margin = margin(r = 12)),
  axis.text = element_text(face = "bold", colour = "black", size = 8),
  legend.position = "bottom",
  panel.spacing = unit(0, "pt"),
  strip.background = element_rect(fill = "black", colour = "goldenrod",
                                  size = 1),
  strip.text = element_text(face = "bold", colour = "white", size = 9)
)

kBoxCountTheme.default <- theme(
  axis.title = element_text(face = "bold", colour = "black", size = 12),
  axis.title.x = element_text(margin = margin(t = 12)),
  axis.text = element_text(face = "bold", colour = "black"),
  axis.text.x = element_text(hjust = 0, vjust = 0.5, angle = 90, size = 7),
  axis.text.y = element_text(size = 10)
)

kSignatureTheme.default <- theme(
  axis.title = element_text(face = "bold", size = 12),
  axis.title.x = element_text(margin = margin(t = 12)),
  axis.title.y = element_text(margin = margin(r = 10)),
  axis.text = element_text(face = "bold", colour = "black"),
  axis.text.y = element_text(family = "sans", size = 7.5),
  strip.background = element_rect(fill = "lavender"),
  strip.text = element_text(face = "bold", size = 9)
)

kContributionTheme.default <- theme(
  axis.title = element_text(face = "bold", size = 12),
  axis.title.x = element_text(margin = margin(t = 12)),
  axis.title.y = element_text(margin = margin(r = 10)),
  axis.text = element_text(face = "bold", colour = "black"),
  axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
  axis.text.y = element_text(size = 10)
)

kCosineTheme.default <- theme(
  axis.title = element_blank(),
  axis.text = element_text(face = "bold", colour = "black", size = 11),
  axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)
)

kMeasureTheme.default <- theme(
  axis.title.x = element_text(face = "bold", size = 12,
                              margin = margin(t = 10)),
  axis.text = element_text(face = "bold", colour = "black"),
  strip.background = element_rect(fill = "lavender"),
  strip.text = element_text(face = "bold", size = 9)
)

kLineGraphTheme.default <- theme(
  axis.title.x = element_text(face = "bold", size = 12,
                              margin = margin(t = 12)),
  axis.title.y = element_blank(),
  axis.text = element_text(face = "bold", colour = "black"),
  axis.text.x = element_text(size = 8, hjust = 1, vjust = 0.5, angle = 90),
  axis.text.y = element_text(size = 10)
)

kLegendTheme.default <- theme(
  legend.background = element_rect(colour = "gray"),
  legend.title = element_text(face = "bold", size = 10.5),
  legend.title.align = 0.5
)

kLegendTheme.hidden <- theme(
  legend.position = "none"
)

# S3 classes carry details of output files and themes to be used by ggplot2
XPlotter <- function(file.out, theme.extra = NULL) {
  if (!is.object(file.out)) {
    stop("Invalid data type of 'file.out'", call. = F)
  }
  
  if (!(class(file.out)[3] %in% kSupportedOutputFormats)) {
    stop("Unrecognised output format", call. = F)
  }
  
  if (is.theme(theme.extra)) {
    theme <- kDefaultTheme + theme.extra
  } else {
    theme <- kDefaultTheme
  }
  
  me <- list(
    file.out = file.out,
    theme = theme
  )
  
  class(me) <- append(class(me), "XPlotter")
  return(me)
}

BarCountPlotter <- function(file.out, use.default.theme = TRUE,
                            theme.extra = NULL) {
  if (use.default.theme) {
    me <- XPlotter(file.out, theme.extra = kBarTheme.default + kBarTheme.num)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "BarCountPlotter")
  return(me)
}

BarPercentagePlotter <- function(file.out, use.default.theme = TRUE,
                                 theme.extra = NULL) {
  if (use.default.theme) {
    me <- XPlotter(file.out, theme.extra = kBarTheme.default + kBarTheme.pct)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "BarPercentagePlotter")
  return(me)
}

StackedBarPlotter <- function(file.out, use.default.theme = TRUE,
                              theme.extra = NULL) {
  if (use.default.theme) {
    me <- XPlotter(file.out, theme.extra = kStackedBarTheme.default +
                     kLegendTheme.default)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "StackedBarPlotter")
  return(me)
}

HeatmapPlotter <- function(file.out, use.default.theme = TRUE,
                           theme.extra = NULL) {
  if (use.default.theme) {
    me <- XPlotter(file.out,
                   theme.extra = kHeatmapTheme.default + kLegendTheme.default)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "HeatmapPlotter")
  return(me)
}

BoxCountPlotter <- function(file.out, use.default.theme = TRUE,
                         theme.extra = NULL) {
  if (use.default.theme) {
    me <- XPlotter(file.out,
                   theme.extra = kBoxCountTheme.default + kLegendTheme.hidden)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "BoxCountPlotter")
  return(me)
}

SignaturePlotter <- function(file.out, use.default.theme = TRUE,
                             theme.extra = NULL) {
  if (use.default.theme) {
    me <- XPlotter(file.out, theme.extra = kSignatureTheme.default)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "SignaturePlotter")
  return(me)
}

ContributionPlotter <- function(file.out, use.default.theme = TRUE,
                                theme.extra = NULL) {
  if (use.default.theme) {
    me <- XPlotter(file.out, theme.extra = kContributionTheme.default +
                     kLegendTheme.default)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "ContributionPlotter")
  return(me)
}

CosinePlotter <- function(file.out, use.default.theme = TRUE,
                          theme.extra = NULL) {
  if (use.default.theme) {
    me <- XPlotter(file.out, theme.extra = kCosineTheme.default)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "CosinePlotter")
  return(me)
}

MeasurePlotter <- function(file.out, use.default.theme = TRUE,
                           theme.extra = NULL) {
  if (use.default.theme) {
    me <- XPlotter(file.out,
                   theme.extra = kMeasureTheme.default + kLegendTheme.hidden)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "MeasurePlotter")
  return(me)
}

LineGraphPlotter <- function(file.out, use.default.theme = TRUE,
                             theme.extra = NULL) {
  if (use.default.theme) {
    me <- XPlotter(file.out,
                   theme.extra = kLineGraphTheme.default + kLegendTheme.default)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "LineGraphPlotter")
  return(me)
}
