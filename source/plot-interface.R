Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
suppressPackageStartupMessages(library(SomaticSignatures))

source("source/transform.R")

PlotSsm6Basics <- function(mut.ctx, plotter, id, geno.name = "",
                           percentage = FALSE) {
  UseMethod("PlotSsm6Basics", mut.ctx)
}

PlotSsm6Spectrum <- function(mut.ctx, plotter, geno.type = "", order = NULL) {
  UseMethod("PlotSsm6Spectrum", mut.ctx)
}

PlotSsm96Heatmap <- function(mut.ctx, plotter, geno.type = "") {
  UseMethod("PlotSsm96Heatmap", mut.ctx)
}

PlotSsmCounts <- function(mut.ctx, plotter, geno.type = "", log10 = TRUE) {
  UseMethod("PlotSsmCounts", mut.ctx)
}

PlotSsmSignatures <- function(ssm.sigs, plotter, geno.type = "") {
  UseMethod("PlotSsmSignatures", ssm.sigs)
}

PlotSsmSigContribution <- function(ssm.sigs, plotter,
                                   geno.type = "", order = NULL) {
  UseMethod("PlotSsmSigContribution", ssm.sigs)
}

PlotCosineSimilarity <- function(ssm.sigs, plotter, geno.type = "") {
  UseMethod("PlotCosineSimilarity", ssm.sigs)
}

PlotMeasures <- function(nmf, plotter, geno.type = "") {
  UseMethod("PlotMeasures", nmf)
}
