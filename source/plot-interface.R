Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
suppressPackageStartupMessages(library(SomaticSignatures))
suppressPackageStartupMessages(library(dendextend))

source("source/transform.R")

PlotSsm6Basics <- function(mut.ctx, plotter, id, geno.name = NULL,
                           percentage = FALSE) {
  if (!is.object(plotter) || !inherits(plotter, "XPlotter")) {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsm6Basics", mut.ctx)
}

PlotSsm6Spectrum <- function(mut.ctx, plotter, geno.type = NULL, order = NULL) {
  if (!is.object(plotter) || !inherits(plotter, "XPlotter")) {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsm6Spectrum", mut.ctx)
}

PlotSsm96Heatmap <- function(mut.ctx, plotter, geno.type = NULL) {
  if (!is.object(plotter) || !inherits(plotter, "XPlotter")) {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsm96Heatmap", mut.ctx)
}

PlotSsmCounts <- function(mut.ctx, plotter, geno.type = NULL, log10 = TRUE) {
  if (!is.object(plotter) || !inherits(plotter, "XPlotter")) {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsmCounts", mut.ctx)
}

PlotSsmSignatures <- function(ssm.sigs, plotter, geno.type = NULL) {
  if (!is.object(plotter) || !inherits(plotter, "XPlotter")) {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsmSignatures", ssm.sigs)
}

PlotSsmSigContribution <- function(ssm.sigs, plotter, geno.type = NULL,
                                   order = NULL) {
  if (!is.object(plotter) || !inherits(plotter, "XPlotter")) {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSsmSigContribution", ssm.sigs)
}

PlotSigPrevalence <- function(prevalence, plotter, n.samples = NULL) {
  if (!is.object(plotter) || !inherits(plotter, "XPlotter")) {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotSigPrevalence", prevalence)
}

PlotCosineSimilarity <- function(ssm.sigs, plotter, geno.type = NULL) {
  if (!is.object(plotter) || !inherits(plotter, "XPlotter")) {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotCosineSimilarity", ssm.sigs)
}

PlotMeasures <- function(nmf, plotter, geno.type = NULL) {
  if (!is.object(plotter) || !inherits(plotter, "XPlotter")) {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotMeasures", nmf)
}

PlotLineComparison <- function(list.stats, plotter, palette = "Set1",
                               title = NULL, log10 = FALSE) {
  if (!is.object(plotter) || !inherits(plotter, "XPlotter")) {
    stop("Invalid 'plotter'", call. = F)
  }
  
  UseMethod("PlotLineComparison", list.stats)
}

PlotDendrogram <- function(fit, file.out, n.clust = NULL,
                           group.frame = TRUE, topic = NULL) {
  if (!is.object(file.out) || !inherits(file.out, "XFileOut")) {
    stop("Invalid 'file.out'", call. = F)
  }
  
  if (!inherits(file.out, "JpgFileOut") && !inherits(file.out, "PdfFileOut")) {
    stop("This type of 'file.out' is currently unsupported", call. = F)
  }
  
  UseMethod("PlotDendrogram", fit)
}
