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
  axis.title.x = element_text(face = "bold", size = 12, margin = margin(15))
)

kBarTheme.pct <- theme(
  axis.title.x = element_blank()
)

kStackedBarTheme.default <- theme(
  axis.title.x = element_text(face = "bold", size = 14, margin = margin(15)),
  axis.title.y = element_blank(),
  axis.text = element_text(face = "bold", colour = "black"),
  axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
  axis.text.y = element_text(size = 11),
  legend.background = element_rect(colour = "gray"),
  legend.title = element_text(face = "bold", size = 10.5),
  legend.title.align = 0.5
)

kHeatmapTheme.default <- theme(
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

# S3 classes carry details of output files and themes to be used by ggplot2
# Four inherited classes: "BarCountPlotter", "BarPercentagePlotter",
# "StackedBarPlotter" and "HeatmapPlotter". Currenly it accepts two types of
# "file.out": "PdfFileOut" and/or "JpgFileOut"
XPlotter <- function(file.out, theme.extra = NA) {
  if (!is.object(file.out)) {
    stop("Invalid data type of 'file.out'", call. = F)
  }
  
  if (!(tail(class(file.out), n = 1) %in% kSupportedOutputFormats)) {
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
                            theme.extra = NA) {
  if (use.default.theme) {
    me <- XPlotter(file.out, theme.extra = kBarTheme.default + kBarTheme.num)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "BarCountPlotter")
  return(me)
}

BarPercentagePlotter <- function(file.out, use.default.theme = TRUE,
                                 theme.extra = NA) {
  if (use.default.theme) {
    me <- XPlotter(file.out, theme.extra = kBarTheme.default + kBarTheme.pct)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "BarPercentagePlotter")
  return(me)
}

StackedBarPlotter <- function(file.out, use.default.theme = TRUE,
                              theme.extra = NA) {
  if (use.default.theme) {
    me <- XPlotter(file.out, theme.extra = kStackedBarTheme.default)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "StackedBarPlotter")
  return(me)
}

HeatmapPlotter <- function(file.out, use.default.theme = TRUE,
                           theme.extra = NA) {
  if (use.default.theme) {
    me <- XPlotter(file.out, theme.extra = kHeatmapTheme.default)
  } else {
    me <- XPlotter(file.out, theme.extra = theme.extra)
  }
  
  class(me) <- append(class(me), "HeatmapPlotter")
  return(me)
}
