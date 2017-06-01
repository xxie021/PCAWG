Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
suppressPackageStartupMessages(library(SomaticSignatures))
suppressPackageStartupMessages(library(dendextend))

source("source/transform.R")

PlotSsm6Basics <- function(mut.ctx, plotter, id, geno.name = "",
                           percentage = FALSE) {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsm6Basics", mut.ctx)
}

PlotSsm6Spectrum <- function(mut.ctx, plotter, geno.type = "", order = NULL) {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsm6Spectrum", mut.ctx)
}

PlotSsm96Heatmap <- function(mut.ctx, plotter, geno.type = "") {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsm96Heatmap", mut.ctx)
}

PlotSsmCounts <- function(mut.ctx, plotter, geno.type = "", log10 = TRUE) {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsmCounts", mut.ctx)
}

PlotSsmSignatures <- function(ssm.sigs, plotter, geno.type = "") {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsmSignatures", ssm.sigs)
}

PlotSsmSigContribution <- function(ssm.sigs, plotter,
                                   geno.type = "", order = NULL) {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsmSigContribution", ssm.sigs)
}

PlotCosineSimilarity <- function(ssm.sigs, plotter, geno.type = "") {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotCosineSimilarity", ssm.sigs)
}

PlotMeasures <- function(nmf, plotter, geno.type = "") {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotMeasures", nmf)
}

PlotDendrogram <- function(fit, file.out, n.clust = NULL,
                           group.frame = TRUE, topic = NULL) {
  if (!is.object(file.out) || class(file.out)[2] != "XFileOut") {
    stop("Invalid 'file.out'", call. = F)
  }
  
  if (!(class(file.out)[3] %in% c("JpgFileOut", "PdfFileOut"))) {
    stop("This type of 'file.out' is currently unsupported", call. = F)
  }
  
  UseMethod("PlotDendrogram", fit)
}
