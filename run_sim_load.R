suppressPackageStartupMessages(library(optparse))

#Constant
kVcfNamePart <- ".consensus.20160830.somatic.indel"

options <- list(
  make_option(c("-l", "--link_file"),
              type = "character",
              help = "Specify a link file stating VCF files per tumour type"),
  make_option(c("-t", "--tumour_type"),
              type = "character",
              help = "Specify tumour type(s) to analyse, separated by comma"),
  make_option(c("-p", "--vcf_path"),
              type = "character",
              default = "",
              help = "Specify a path to access VCF files"),
  make_option(c("-d", "--directory"),
              type = "character",
              default = "sim_data/",
              help = paste("Specify a directory to store SIM count data",
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
  sim.set <- mapply(function(fname, id, name) {
    tryCatch({
      cat("Info: Loading VCF file \"", basename(fname), "\" ...\n", sep = "")
      raw.vcf <- readVcf(fname, "hg19", ScanVcfParam(info = NA, samples = NA))
      sim <- SimSample(raw.vcf, id = id, name = name)
      cat("Info: Loading completed\n")
      return(sim)
    }, error = function(e) {
      stop("Failed to load this VCF file. Stop", call. = F)
    })
  },
  fname = paste0(opt$vcf_path, f$tumor_wgs_aliquot_id, kVcfNamePart, ".vcf.gz"),
  id = f$t_sample_index,
  name = f$tumor_wgs_aliquot_id,
  SIMPLIFY = F)
  
  if (length(sim.set) == 1) {
    mut.len <- sim.set[[1]]$mut.ctx
  } else {
    mut.len <- MergeMotifMatrixFromSamples(sim.set)
    mut.len <- mut.len[, order(as.integer(colnames(mut.len)))]
  }
  geno.type <- tolower(unique(f$histology_abbreviation))
  cat("Info: Preparing output directory ...\n")
  data.length.fo <- RDataFileOut(paste0(geno.type, ".indel.length.counts"),
                                 out.path = opt$directory)
  Save2RData(mut.len, data.length.fo,
             obj.name = paste0(stringr::str_replace_all(geno.type, "-", "."),
                               ".mut.len"))
})
