Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
suppressPackageStartupMessages(library(SomaticSignatures))

# S3 classes store variant info of a sample
VariantSample <- function(vcf, ref, id = "", name = "") {
  if (class(vcf) != "CollapsedVCF") {
    stop("Invalid VCF data found", call. = F)
  }
  
  if (class(ref) != "FaFile") {
    stop("Invalid reference sequence found", call. = F)
  }
  
  #tryCatch({
  #  countFa(ref)
  #}, error = function(e) {
  #  stop("Invalid reference sequence found", call. = F)
  #})
  
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
    ref = ref,
    id = id,
    name = name
  )
  
  class(me) <- append(class(me), "VariantSample")
  return(me)
}

# Somatic substitution mutation
SsmSample <- function(vcf, ref, id = "", name = "") {
  class.name <- "SsmSample"
  
  cat("Info: Constructing class \"", class.name, "\" ...\n", sep = "")
  me <- VariantSample(vcf, ref, id, name)
  
  idx <- which(isSNV(me$vcf))
  seq.ctx <- mutationContext(as(vcf[idx], "VRanges"), ref)
  sampleNames(seq.ctx) <- GenSampleName(me)
  
  me$idx <- idx
  me$mut.ctx <- motifMatrix(seq.ctx, normalize = F)
  me$mut.total <- sum(me$mut.ctx)
  
  class(me) <- append(class(me), class.name)
  cat("Info: Class \"", class.name, "\" constructed\n", sep = "")
  
  return(me)
}

# Somatic indel mutation
# SimSample <- function(vcf, ref, id = "", name = "") {
#   class.name <- "SimSample"
#   
#   cat("Info: Constructing class \"", class.name, "\" ...\n", sep = "")
#   me <- VariantSample(vcf, ref, id, name)
#   idx.ins <- which(isInsertion(me$vcf))
#   idx.del <- which(isDeletion(me$vcf))
#   
#   me$idx.ins <- idx.ins
#   me$idx.del <- idx.del
#   
#   class(me) <- append(class(me), class.name)
#   cat("Info: Class \"", class.name, "\" constructed\n", sep = "")
#   
#   return(me)
# }

GenSampleName <- function(sample, full = FALSE) {
  UseMethod("GenSampleName", sample)
}

GenSampleName.VariantSample <- function(sample, full = FALSE) {
  if (full && trimws(sample$id) != "" && trimws(sample$name) != "") {
    return(paste0(sample$id, ".", sample$name))
  }
  
  return(ifelse(trimws(sample$id) == "", sample$name, sample$id))
}
