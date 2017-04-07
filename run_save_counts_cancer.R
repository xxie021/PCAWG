args <- commandArgs(trailingOnly = T)
#args <- c("", "vcf.gz", "refseq.fasta", "pdf", "Lung", "Ovary")

# Constants
kMinArgLength <- 3

if (length(args) < kMinArgLength) {
  stop("Missing arguments", call. = F)
}

tryCatch({
  # Mandatory: load RData for the relationship between samples and tumour types
  load("tsample2ttype.RData")
}, error = function(e) {
  stop("Failed to load the \"tsample2ttype.RData\" file", call. = F)
})

if (length(args) == kMinArgLength) {
  files <- tsample2ttype
} else {
  files <- tsample2ttype[tsample2ttype$histology_abbreviation %in%
                           args[kMinArgLength + 1:length(args)], ]
}

if (length(files$tumor_wgs_aliquot_id) == 0) {
  stop("No associated VCF files found", call. = F)
}

if (!exists("VariantSet", mode = "function")) {
  source("variant_set.R")
}

files.per.tumour <- split(files, f = files$histology_abbreviation)
bin <- lapply(files.per.tumour, function(f) {
  singlessm.set <- mapply(function(fname, id, name) {
    tryCatch({
      cat("Loading VCF file: \"", basename(fname), "\" ...\n", sep = "")
      raw.vcf <- readVcf(fname, "hg19", ScanVcfParam(info = NA, samples = NA))
      singlessm <- SingleSSM(raw.vcf, FaFile(args[3]), id = id, name = name)
      cat("Loading completed\n")
      return(singlessm)
    }, error = function(e) {
      stop("Failed to load this VCF file. Stop", call. = F)
    })
  },
  fname = paste0(args[1], f$tumor_wgs_aliquot_id, ".", args[2]),
  id = f$t_sample_index,
  name = f$tumor_wgs_aliquot_id,
  SIMPLIFY = F)
  
  SaveMutationCounts(SsmSampleSet(singlessm.set, f$histology_abbreviation[1]))
})
