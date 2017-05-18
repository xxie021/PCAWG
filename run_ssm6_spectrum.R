suppressPackageStartupMessages(library(optparse))

source("source/models/file-out.R")

options <- list(
  make_option(c("-m", "--motif_matrix"),
              type = "character",
              help = "Specify a motif matrix to process"),
  make_option(c("-d", "--parent_dir"),
              type = "character",
              help = "Specify a directory to store spectrum results"),
  make_option(c("-s", "--sort_by"),
              type = "character",
              help = paste("Specify one of the six base substitutions to",
                           "order result")),
  make_option(c("--jpg"),
              action = "store_true",
              default = FALSE,
              help = "Indicate if storing spectrum plots in JPG format")
)

opt <- tryCatch({
  parse_args(OptionParser(option_list = options))
}, warning = function(w) {
}, error = function(e) {
  stop("Incorrect usage. Please see --help", call. = F)
})

start.time <- Sys.time()
cat("Info: Starting script at ", strftime(start.time), " ...\n", sep = "")

if (is.null(opt$motif_matrix)) {
  stop("MOTIF_MATRIX must be provided", call. = F)
} else {
  geno.type <- strsplit(basename(opt$motif_matrix), "\\.")[[1]][1]
}

if (is.null(opt$parent_dir)) {
  opt$parent_dir <- paste0(strftime(start.time, format = "%Y%m%d%H%M%S"), "/")
} else {
  opt$parent_dir <- trimws(opt$parent_dir)
  if (opt$parent_dir != "" && !endsWith(opt$parent_dir, "/") &&
      !endsWith(opt$parent_dir, "\\")) {
    opt$parent_dir <- paste0(opt$parent_dir, "/")
  }
}

cat("Info: Preparing output directory ...\n")
data.spectrum.fo <- TsvFileOut(paste0(geno.type, ".ssm6.spectrum"),
                               out.path = paste0(opt$parent_dir, "data"))
if (opt$jpg) {
  plot.stackedbar.fo <- JpgFileOut(paste0(geno.type, ".ssm6.spectrum"),
                                   out.path = paste0(opt$parent_dir, "plots"))
} else {
  plot.stackedbar.fo <- PdfFileOut(paste0(geno.type, ".ssm6.spectrum"),
                                   out.path = paste0(opt$parent_dir, "plots"))
}

cat("Info: Loading data \"", basename(opt$motif_matrix), "\" ...\n", sep = "")
mut.ctx <- tryCatch({
  get(load(opt$motif_matrix))
}, warning = function(w) {
  stop("MOTIF_MATRIX RData file can't be loaded", call. = F)
}, error = function(e) {
  stop("MOTIF_MATRIX RData file can't be loaded", call. = F)
})
if (!is.matrix(mut.ctx) || nrow(mut.ctx) != 96) {
  stop("Invalid motif matrix. Must be a matrix with 96 rows (mutation types)",
       call. = F)
}

cat("Info: Loading libraries and scripts ...\n")
source("source/models/plotter.R")
source("source/plot-core.R")
source("source/summary.R")

if (geno.type == "wg") {
  geno.type = "Whole Genomes"
}
stackedbar.plotter <- StackedBarPlotter(plot.stackedbar.fo)
PlotSsm6Spectrum(mut.ctx, stackedbar.plotter, geno.type = geno.type,
                 order = opt$sort_by)
Summary(mut.ctx, "base.summary", file.out = data.spectrum.fo)
