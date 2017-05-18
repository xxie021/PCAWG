suppressPackageStartupMessages(library(optparse))

source("source/models/file-out.R")

options <- list(
  make_option(c("-m", "--motif_matrix"),
              type = "character",
              help = "Specify a motif matrix to process"),
  make_option(c("-n", "--normalise"),
              action = "store_true",
              default = FALSE,
              help = "Indicate if the motif matrix is normalised to run NMF"),
  make_option(c("-i", "--iteration"),
              type = "integer",
              default = 500,
              help = paste("Specify a positive number of iterations",
                           "to run NMF [default %default]")),
  make_option(c("-k", "--n_signatures"),
              type = "integer",
              default = 0,
              help = paste("Specify a nonnegative number of signatures",
                           "[default %default]")),
  make_option(c("-s", "--sort_contribution_by"),
              type = "character",
              help = "Specify a signature for contribution sorting"),
  make_option(c("-d", "--parent_dir"),
              type = "character",
              help = "Specify a parent directory to store signature results"),
  make_option(c("--cosine_low"),
              type = "double",
              default = 0.55,
              help = paste("Specify a threshold (0 < l < 1) of cosine",
                           "similarity to indicate two signatures", 
                           "are different [default %default]")),
  make_option(c("--cosine_high"),
              type = "double",
              default = 0.85,
              help = paste("Specify a threshold (0 < h < 1) of cosine",
                           "similarity to indicate two signatures are", 
                           "redundant [default %default]")),
  make_option(c("--jpg"),
              action = "store_true",
              default = FALSE,
              help = "Indicate if storing plots in JPG format"),
  
  make_option(c("--log10_count"),
              action = "store_true",
              default = FALSE,
              help = paste("Indicate if plotting mutation type counts",
                           "using Log 10 scale"))
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

if (is.null(opt$iteration) || opt$iteration < 1) {
  stop("ITERATION must be a positive number", call. = F)
}

if (is.null(opt$n_signatures) || opt$n_signatures < 0) {
  stop("N_SIGNATURES must be a nonnegative number", call. = F)
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

if (is.null(opt$cosine_low) || is.null(opt$cosine_high) ||
    opt$cosine_low <= 0 || opt$cosine_high >= 1 ||
    opt$cosine_low > opt$cosine_high) {
  stop("COSINE threahold must be a fraction satisfying (0 < low <= high < 1)",
       call. = F)
}

cat("Info: Preparing output directories ...\n")
data.nmf.fo <- RDataFileOut(paste0(geno.type, ".ssm.nmf"),
                            out.path = paste0(opt$parent_dir, "data"))
data.sig.fo <- RDataFileOut(paste0(geno.type, ".ssm.sigs"),
                            out.path = paste0(opt$parent_dir, "data"))
if (opt$jpg) {
  plot.count.fo <- JpgFileOut(paste0(geno.type, ".ssm.counts"),
                              out.path = paste0(opt$parent_dir, "plots"))
  plot.sig.fo <- JpgFileOut(paste0(geno.type, ".ssm.sigs"),
                            out.path = paste0(opt$parent_dir, "plots"))
  plot.contri.fo <- JpgFileOut(paste0(geno.type, ".ssm.contri"),
                               out.path = paste0(opt$parent_dir, "plots"))
  plot.cos.fo <- JpgFileOut(paste0(geno.type, ".ssm.cos"),
                            out.path = paste0(opt$parent_dir, "plots"))
  if (opt$n_signatures < 2) {
    plot.meas.fo <- JpgFileOut(paste0(geno.type, ".ssm.meas"),
                               out.path = paste0(opt$parent_dir, "plots"))
  }
} else {
  plot.count.fo <- PdfFileOut(paste0(geno.type, ".ssm.counts"),
                              out.path = paste0(opt$parent_dir, "plots"))
  plot.sig.fo <- PdfFileOut(paste0(geno.type, ".ssm.sigs"),
                            out.path = paste0(opt$parent_dir, "plots"))
  plot.contri.fo <- PdfFileOut(paste0(geno.type, ".ssm.contri"),
                               out.path = paste0(opt$parent_dir, "plots"))
  plot.cos.fo <- PdfFileOut(paste0(geno.type, ".ssm.cos"),
                            out.path = paste0(opt$parent_dir, "plots"))
  if (opt$n_signatures < 2) {
    plot.meas.fo <- PdfFileOut(paste0(geno.type, ".ssm.meas"),
                               out.path = paste0(opt$parent_dir, "plots"))
  }
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

if (ncol(mut.ctx) < 10) {
  stop("Too few samples/genomes", call. = F)
}

cat("Info: Loading libraries and scripts ...\n")
source("source/models/plotter.R")
source("source/ssm-core.R")
source("source/plot-core.R")
source("source/utils.R")

cat("Info: Removing minor mutation types ...\n")
mut.ctx.reduced <- RemoveMinorMutationTypes(mut.ctx, threshold = 0.01,
                                            normalise = opt$normalise)

if (opt$n_signatures > min(nrow(mut.ctx.reduced), ncol(mut.ctx.reduced)) - 1) {
  stop("Too large N_SIGNATURES specified", call. = F)
}

cat("Info: Bootstrapping ...\n")
mut.ctx.boot <- McBootstrap(mut.ctx.reduced)

if (opt$normalise) {
  cat("Info: Normalising data ...\n")
  mut.ctx.boot <- sweep(mut.ctx.boot, 2, colSums(mut.ctx.boot), `/`)
}

if (opt$n_signatures < 2) {
  # Number of signatures is not provided
  cat("Info: Estimating 'k' range ...\n")
  max.k <- EstMaxNumOfSignatures(ncol(mut.ctx.boot))
  cat("Info: 'k' in [2, ", max.k, "] will be tested\n", sep = "")

  if (max.k > 2) {
    cat("Info: Running NMF ...\n")
    nmf <- nmf(mut.ctx.boot, 2:max.k, "brunet", seed = 42, nrun = opt$iteration)
    Save2RData(nmf, data.nmf.fo,
               obj.name = paste0(stringr::str_replace_all(tolower(geno.type),
                                                          "-", "."), ".nmf"))
    cat("Info: NMF completed\n")

    sigs <- ConstSignatures(nmf, mut.ctx.boot, cosine.low = opt$cosine_low,
                            cosine.high = opt$cosine_high)
  } else {
    cat("Info: Estimated number of signatures: 2\n")
    sigs <- identifySignatures(mut.ctx.boot, 2, nmfDecomposition, seed = 42)
  }
} else {
  # Number of signatures is provided
  cat("Info: Input number of signatures: ", opt$n_signatures, "\n", sep = "")
  sigs <- identifySignatures(mut.ctx.boot, opt$n_signatures, nmfDecomposition,
                             seed = 42)
}
Save2RData(sigs, data.sig.fo,
           obj.name = paste0(stringr::str_replace_all(tolower(geno.type),
                                                      "-", "."), ".sigs"))

if (geno.type == "wg") {
  geno.type = "Whole Genomes"
}

count.plotter <- BoxCountPlotter(plot.count.fo)
PlotSsmCounts(mut.ctx, count.plotter, geno.type = geno.type, opt$log10_count)

sig.plotter <- SignaturePlotter(plot.sig.fo)
PlotSsmSignatures(sigs, sig.plotter, geno.type = geno.type)

contri.plotter <- ContributionPlotter(plot.contri.fo)
PlotSsmSigContribution(sigs, contri.plotter, geno.type = geno.type,
                       order = opt$sort_contribution_by)

cos.plotter <- CosinePlotter(plot.cos.fo)
PlotCosineSimilarity(sigs, cos.plotter, geno.type = geno.type)

if (opt$n_signatures < 2 && max.k > 2) {
  meas.plotter <- MeasurePlotter(plot.meas.fo)
  PlotMeasures(nmf, meas.plotter, geno.type = geno.type)
}
