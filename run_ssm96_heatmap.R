suppressPackageStartupMessages(library(optparse))

options <- list(
  make_option(c("-m", "--motif_matrix"),
              type = "character",
              help = "Specify a motif matrix to process"),
  make_option(c("-d", "--directory"),
              type = "character",
              help = "Specify a directory to store heatmap"),
  make_option(c("--jpg"),
              action = "store_true",
              default = FALSE,
              help = "Indicate if storing heatmap plot in JPG format")
)

opt <- tryCatch({
  parse_args(OptionParser(option_list = options))
}, warning = function(w) {
}, error = function(e) {
  stop("Incorrect usage. Please see --help", call. = F)
})

start.time <- Sys.time()
cat("Info: Starting script at ", strftime(start.time), " ...\n", sep = "")

cat("Info: Loading libraries and scripts ...\n")
source("source/models/file-out.R")
source("source/models/plotter.R")
source("source/plot-core.R")

if (is.null(opt$motif_matrix)) {
  stop("MOTIF_MATRIX must be provided", call. = F)
} else {
  geno.type <- strsplit(basename(opt$motif_matrix), "\\.")[[1]][1]
}

if (is.null(opt$directory)) {
  opt$directory <- paste0(strftime(start.time, format = "%Y%m%d%H%M%S"), "/")
} else {
  opt$directory <- trimws(opt$directory)
  if (opt$directory != "" && !endsWith(opt$directory, "/") &&
      !endsWith(opt$directory, "\\")) {
    opt$directory <- paste0(opt$directory, "/")
  }
}

cat("Info: Preparing output directory ...\n")
if (opt$jpg) {
  heatmap.fo <- JpgFileOut(paste0(geno.type, ".ssm96.heatmap"),
                           out.path = opt$directory)
} else {
  heatmap.fo <- PdfFileOut(paste0(geno.type, ".ssm96.heatmap"),
                           out.path = opt$directory)
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

heatmap.plotter <- HeatmapPlotter(heatmap.fo)

if (geno.type == "wg") {
  geno.type = "Whole Genomes"
}
PlotSsm96Heatmap(mut.ctx, heatmap.plotter, geno.type = geno.type)
