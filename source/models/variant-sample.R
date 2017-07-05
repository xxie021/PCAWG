Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
suppressPackageStartupMessages(library(SomaticSignatures))

# S3 classes store variant info of a sample
VariantSample <- function(vcf, id = "", name = "") {
  if (class(vcf) != "CollapsedVCF") {
    stop("Invalid VCF data found", call. = F)
  }
  
  if (!is.character(id)) {
    cat("Warn: 'id' is not type of character. Converted\n")
    id <- as.character(id)
  }
  
  if (!is.character(name)) {
    cat("Warn: 'name' is not type of character. Converted\n")
    name <- as.character(name)
  }
  
  if (trimws(id) == "" && trimws(name) == "") {
    stop("'id' and 'name' can't be both empty", call. = F)
  }
  
  me <- list(
    vcf = vcf,
    id = id,
    name = name
  )
  
  class(me) <- append(class(me), "VariantSample")
  return(me)
}

# Somatic substitution mutation
SsmSample <- function(vcf, ref, id = "", name = "") {
  if (class(ref) != "FaFile") {
    stop("Invalid reference sequence found", call. = F)
  }
  
  #tryCatch({
  #  countFa(ref)
  #}, error = function(e) {
  #  stop("Invalid reference sequence found", call. = F)
  #})
  
  class.name <- "SsmSample"
  
  cat("Info: Constructing class \"", class.name, "\" ...\n", sep = "")
  me <- VariantSample(vcf, id, name)
  
  idx <- which(isSNV(me$vcf))
  seq.ctx <- mutationContext(as(vcf[idx], "VRanges"), ref)
  sampleNames(seq.ctx) <- GenSampleName(me)
  
  me$ref <- ref
  me$idx <- idx
  me$mut.ctx <- motifMatrix(seq.ctx, normalize = F)
  me$mut.total <- sum(me$mut.ctx)
  
  class(me) <- append(class(me), class.name)
  cat("Info: Class \"", class.name, "\" constructed\n", sep = "")
  
  return(me)
}

# Somatic indel mutation
SimSample <- function(vcf, id = "", name = "") {
  class.name <- "SimSample"

  cat("Info: Constructing class \"", class.name, "\" ...\n", sep = "")
  me <- VariantSample(vcf, id, name)
  idx.ins <- which(isInsertion(me$vcf))
  idx.del <- which(isDeletion(me$vcf))

  me$idx.ins <- idx.ins
  me$idx.del <- idx.del
  
  len.all <- sapply(which(isIndel(me$vcf)), function(idx) {
    len <- width(alt(me$vcf)[[idx]]) - width(ref(me$vcf)[idx])
    if (len == 1) {
      ins.base <- strsplit(as.character(alt(me$vcf)[[idx]]), split = "")[[1]][2]
      len <- ifelse(ins.base %in% c("C", "G"), 0.67, 0.84)
    } else if (len == -1) {
      del.base <- strsplit(as.character(ref(me$vcf)[idx]), split = "")[[1]][2]
      len <- ifelse(del.base %in% c("C", "G"), -0.67, -0.84)
    }
    return(len)
  })
  
  me$mut.len <- matrix(c(length(which(len.all == 0.67)),
                         length(which(len.all == 0.84)),
                         length(which(len.all > 1 & len.all <= 5)),
                         length(which(len.all > 5)),
                         length(which(len.all == -0.67)),
                         length(which(len.all == -0.84)),
                         length(which(len.all < -1 & len.all >= -5)),
                         length(which(len.all < -5))), ncol = 1,
                       dimnames = list(c("ins.C", "ins.T", "ins.2to5",
                                         "ins.5p", "del.C", "del.T", "del.2to5",
                                         "del.5p"), GenSampleName(me)))

  class(me) <- append(class(me), class.name)
  cat("Info: Class \"", class.name, "\" constructed\n", sep = "")

  return(me)
}

GenSampleName <- function(sample, full = FALSE) {
  UseMethod("GenSampleName", sample)
}

GenSampleName.VariantSample <- function(sample, full = FALSE) {
  if (full && trimws(sample$id) != "" && trimws(sample$name) != "") {
    return(paste0(sample$id, ".", sample$name))
  }
  
  return(ifelse(trimws(sample$id) == "", sample$name, sample$id))
}
