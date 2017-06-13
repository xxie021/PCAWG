suppressPackageStartupMessages(library(optparse))

source("source/models/file-out.R")

options <- list(
  make_option(c("-s", "--sample_matrix"),
              type = "character",
              help = paste("Specify a sample matrix or a directory of",
                           "sample matrices to process")),
  make_option(c("--all"),
              action = "store_true",
              default = FALSE,
              help = paste("Indicate if '-s' specifies a directory and",
                           "load all files in that directory")),
  make_option(c("-c", "--consensus_signatures"),
              type = "character",
              help = paste("Specify a consensus signatures matrix to",
                           "perform reconstruction. If not specified,",
                           "COSMIC published signatures are to be used")),
  make_option(c("-t", "--threshold"),
              type = "double",
              default = 0.25,
              help = paste("Specify a threshold to filter out insignificant",
                           "mutational signatures [default %default]")),
  make_option(c("-d", "--parent_dir"),
              type = "character",
              help = "Specify a parent directory to store prevalence results"),
  make_option(c("--tsv"),
              action = "store_true",
              default = FALSE,
              help = "Indicate if storing data in TSV format")
)

opt <- tryCatch({
  parse_args(OptionParser(option_list = options))
}, warning = function(w) {
}, error = function(e) {
  stop("Incorrect usage. Please see --help", call. = F)
})

start.time <- Sys.time()
cat("Info: Starting script at ", strftime(start.time), " ...\n", sep = "")

if (is.null(opt$sample_matrix) || trimws(opt$sample_matrix) == "") {
  stop("SAMPLE_MATRIX must be provided", call. = F)
} else if (opt$all) {
  sample.files <- list.files(path = opt$sample_matrix, pattern = "[.]RData$",
                             full.names = T)
} else {
  sample.files <- opt$sample_matrix
}

if (length(sample.files) == 0) {
  stop("No sample matrix is found", call.F = F)
}

if (opt$threshold <= 0 || opt$threshold > 1) {
  stop("THRESHOLD must be a positive fraction not greater than 1", call. = F)
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

cat("Info: Preparing output directories ...\n")
sigs.type <- ifelse(is.null(opt$consensus_signatures), "cosmic",
                    tools::file_path_sans_ext(
                      basename(opt$consensus_signatures)))
if (opt$tsv) {
  data.summary.fo <- TsvFileOut(paste0(sigs.type, ".prevalence"),
                                out.path = paste0(opt$parent_dir, "summary"))
} else {
  data.summary.fo <- RDataFileOut(paste0(sigs.type, ".prevalence"),
                                  out.path = paste0(opt$parent_dir, "summary"))
}

sigs.consensus <- NULL
if (!is.null(opt$consensus_signatures) &&
    trimws(opt$consensus_signatures) != "") {
  cat("Info: Loading consensus signatures \"",
      basename(opt$consensus_signatures), "\" ...\n", sep = "")
  sigs.consensus <- tryCatch({
    get(load(opt$consensus_signatures))
  }, warning = function(w) {
    stop("CONSENSUS_SIGNATURES RData file can't be loaded", call. = F)
  }, error = function(e) {
    stop("CONSENSUS_SIGNATURES RData file can't be loaded", call. = F)
  })
  if (!is.matrix(sigs.consensus) || nrow(sigs.consensus) != 96) {
    stop("Invalid consensus signatures matrix. Must be a matrix with 96 rows",
         call. = F)
  }
}

cat("Info: Loading libraries and scripts ...\n")
source("source/ssm-core.R")
source("source/utils.R")

sigs.prevalence <- do.call(rbind, lapply(sample.files, function(file) {
  cat("Info: Loading data \"", basename(file), "\" ...\n", sep = "")
  mt.samples <- tryCatch({
    get(load(file))
  }, warning = function(w) {
    cat("Warn: SAMPLE_MATRIX RData file can't be loaded. Skipped\n")
    return(NULL)
  }, error = function(e) {
    cat("Warn: SAMPLE_MATRIX RData file can't be loaded. Skipped\n")
    return(NULL)
  })
  if (!is.null(mt.samples) &&
      (!is.matrix(mt.samples) || nrow(mt.samples) != 96)) {
    cat("Warn: Invalid sample matrix (less than 96 rows). Skipped\n")
    return(NULL)
  }
  
  if (!is.null(mt.samples)) {
    tumour.type <- strsplit(basename((file)), split = "\\.")[[1]][1]
    res <- SummarizeSsm96SignaturePrevalence(mt.samples, tumour.type,
                                             mt.sigs = sigs.consensus,
                                             threshold = opt$threshold)
    
    if (opt$tsv) {
      data.raw.fo <- TsvFileOut(paste0(tumour.type, ".reconst"),
                                out.path = paste0(opt$parent_dir, "raw"))
      cat("Info: Saving file \"", basename(data.raw.fo$fullname), "\" ...\n",
          sep = "")
      write.table(res$raw, data.raw.fo$fullname, quote = F, sep = "\t",
                  col.names = NA)
      cat("Info: File \"", basename(data.raw.fo$fullname), "\" saved\n",
          sep = "")
    } else {
      data.raw.fo <- RDataFileOut(paste0(tumour.type, ".reconst"),
                                  out.path = paste0(opt$parent_dir, "raw"))
      Save2RData(res$raw, data.raw.fo,
                 obj.name = paste0(stringr::str_replace_all(tumour.type, "-", "."),
                                   ".reconst"))
    }
    
    return(res$summary)
  }
  
  return(NULL)
}))

if (!is.null(sigs.prevalence)) {
  if (opt$tsv) {
    cat("Info: Saving file \"", basename(data.summary.fo$fullname), "\" ...\n",
        sep = "")
    write.table(sigs.prevalence, data.summary.fo$fullname, quote = F,
                sep = "\t", col.names = NA)
    cat("Info: File \"", basename(data.summary.fo$fullname), "\" saved\n",
        sep = "")
  } else {
    Save2RData(sigs.prevalence, data.summary.fo,
               obj.name = paste0(stringr::str_replace_all(sigs.type, "-", "."),
                                 ".prevalence"))
  }
} else {
  cat("Warn: All sample matrices are invalid and skipped\n")
}
