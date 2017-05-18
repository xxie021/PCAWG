suppressPackageStartupMessages(library(optparse))

#Constant
kVcfNamePart <- ".consensus.20160830.somatic.snv_mnv"

options <- list(
  make_option(c("-l", "--link_file"),
              type = "character",
              help = "Specify a link file stating VCF files per tumour type"),
  make_option(c("-t", "--tumour_type"),
              type = "character",
              help = "Specify tumour type(s) to analyse, separated by comma"),
  make_option(c("-r", "--refseq"),
              type = "character",
              help = "Specify a reference sequence file"),
  make_option(c("-p", "--vcf_path"),
              type = "character",
              default = "",
              help = "Specify a path to access VCF files"),
  make_option(c("-d", "--directory"),
              type = "character",
              default = "data/",
              help = paste("Specify a directory to store SSM96 count data",
                           "[default %default]"))
)

opt <- tryCatch({
  parse_args(OptionParser(option_list = options))
}, warning = function(w) {
}, error = function(e) {
  stop("Incorrect usage. Please see --help", call. = F)
})

start.time <- Sys.time()
cat("Info: Starting script at ", strftime(start.time), " ...\n", sep = "")

if (is.null(opt$link_file)) {
  stop("LINK_FILE must be provided", call. = F)
}

if (is.null(opt$refseq)) {
  stop("REFSEQ must be provided", call. = F)
}

opt$vcf_path <- trimws(opt$vcf_path)
if (opt$vcf_path != "" && !endsWith(opt$vcf_path, "/") &&
    !endsWith(opt$vcf_path, "\\")) {
  opt$vcf_path <- paste0(opt$vcf_path, "/")
}

opt$directory <- trimws(opt$directory)
if (opt$directory != "" && !endsWith(opt$directory, "/") &&
    !endsWith(opt$directory, "\\")) {
  opt$directory <- paste0(opt$directory, "/")
}

cat("Info: Loading link file \"", basename(opt$link_file), "\" ...\n", sep = "")
link.file <- tryCatch({
  get(load(opt$link_file))
}, warning = function(w) {
  stop("LINK_FILE RData file can't be loaded", call. = F)
}, error = function(e) {
  stop("LINK_FILE RData file can't be loaded", call. = F)
})

if (is.character(opt$tumour_type)) {
  cat("Info: Searching VCF files with the specified type(s) ...\n")
  opt$tumour_type <- trimws(unlist(strsplit(opt$tumour_type, split = ",")))
  link.file <- link.file[
    tolower(link.file$histology_abbreviation) %in% tolower(opt$tumour_type),
  ]
}

if (length(link.file$tumor_wgs_aliquot_id) == 0) {
  stop("No relevant VCF files are found", call. = F)
}

cat("Info: Loading libraries and scripts ...\n")
source("source/models/variant-sample.R")
source("source/models/file-out.R")
source("source/utils.R")

files.per.tumour <- split(link.file, f = link.file$histology_abbreviation)
bin <- lapply(files.per.tumour, function(f) {
  cat("Info: Loading tumour type \"", unique(f$histology_abbreviation), "\", ",
      "total #sample(s): ", length(f$histology_abbreviation), "\n",
      sep = "")
  ssm.set <- mapply(function(fname, id, name) {
    tryCatch({
      cat("Info: Loading VCF file \"", basename(fname), "\" ...\n", sep = "")
      raw.vcf <- readVcf(fname, "hg19", ScanVcfParam(info = NA, samples = NA))
      ssm <- SsmSample(raw.vcf, FaFile(opt$refseq), id = id, name = name)
      cat("Info: Loading completed\n")
      return(ssm)
    }, error = function(e) {
      stop("Failed to load this VCF file. Stop", call. = F)
    })
  },
  fname = paste0(opt$vcf_path, f$tumor_wgs_aliquot_id, kVcfNamePart, ".vcf.gz"),
  id = f$t_sample_index,
  name = f$tumor_wgs_aliquot_id,
  SIMPLIFY = F)
  
  if (length(ssm.set) == 1) {
    mut.ctx <- ssm.set[[1]]$mut.ctx
  } else {
    mut.ctx <- MergeMotifMatrixFromSample(ssm.set)
    mut.ctx <- mut.ctx[, order(colnames(mut.ctx))]
  }
  geno.type <- tolower(unique(f$histology_abbreviation))
  cat("Info: Preparing output directory ...\n")
  data.count.fo <- RDataFileOut(paste0(geno.type, ".mut.ctx3.counts"),
                                out.path = opt$directory)
  Save2RData(mut.ctx, data.count.fo,
             obj.name = paste0(stringr::str_replace_all(geno.type, "-", "."),
                               ".mut.ctx"))
})
