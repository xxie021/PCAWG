# Merges multiple motif matrices. If they are to be merged by tumour types,
# i.e. each set of motif matrix represents a specific type of tumour, all
# sample-based counts will be added up prior to merging.
# @param  mut.ctx.set   a list of motif matrices
# @param  by.tumour     the flag indicating if merging is to be done by
#                       tumour types
# @param  tumour.types  an array of tumour types having the same length of 
#                       {@code mut.ctx.set}; used only when {@code by.tumour}
#                       is set to TRUE
# @return               the merged motif matrix
MergeMotifMatrix <- function(mut.ctx.set, by.tumour = FALSE,
                             tumour.types = "") {
  if (!is.list(mut.ctx.set) || length(mut.ctx.set) < 1) {
    stop("Invalid or empty motif matrix set", call. = F)
  }
  
  # Don't know why sometimes sapply returns a list.
  # Call unlist to firmly get a vector
  if (any(unlist(sapply(mut.ctx.set, nrow)) != nrow(mut.ctx.set[[1]]))) {
    stop("Motif matrices in the set have different number of rows", call. = F)
  }
  
  if (by.tumour) {
    if (length(mut.ctx.set) != length(tumour.types)) {
      stop("# motif matrices is not equal to # tumour types", call. = F)
    }
    merged <- do.call(cbind, lapply(mut.ctx.set, function(mut.ctx) {
      rowSums(mut.ctx)
    }))
    colnames(merged) <- tumour.types
  } else {
    if (length(mut.ctx.set) < 2) {
      stop("Single motif matrix is found. Nothing to merge", call. = F)
    }
    merged <- do.call(cbind, mut.ctx.set)
  }
  
  return(merged)
}

# Merges multiple SSM or SIM motif matrices from the given variant sample set.
# @param  sample.set  a list of the {@code SsmSample} or
#                     {@code SimSample} objects
# @return             the merged motif matrix
MergeMotifMatrixFromSamples <- function(sample.set) {
  if (!is.list(sample.set) || length(sample.set) < 2) {
    stop("Invalid or too few VariantSample set", call. = F)
  }
  
  if (!all(sapply(sample.set, inherits, "SsmSample")) &&
      !all(sapply(sample.set, inherits, "SimSample"))) {
    stop("Invalid or different types of VariantSample are found in the set",
         call. = F)
  }
  
  return(MergeMotifMatrix(lapply(sample.set, function(sample) {
    if (inherits(sample, "SsmSample")) {
      return(sample$mut.ctx)
    }
    
    return(sample$mut.len)
  })))
}

# Indicates the tumour type(s) that each consensus mutational signature
# comes from.
# @param  clusters    a named vector indicating the cluster indexes for each
#                     tumour type (names), usually the return of the
#                     {@code cutree} function
# @param  sig.prefix  the string name prefix for consensus mutational signatures
# @return             the matrix showing the tumour type(s) that each consensus
#                     mutational signature comes from
LinkTumourType2Signature <- function(clusters, sig.prefix = "S") {
  indv.names <- names(clusters)
  if (is.null(indv.names)) {
    stop("Invalid 'clusters'. Must be the return of the 'cutree' function",
         call. = F)
  }
  
  if (!is.character(sig.prefix) || length(sig.prefix) != 1 ||
      trimws(sig.prefix) == "") {
    cat("Warn: 'sig.prefix' is incorrect. Use default\n")
    sig.prefix <- "S"
  }
  sigs.consensus <- paste0(sig.prefix, ".", sort(unique(clusters)))
  
  tumour.types <- unique(sapply(indv.names, function(name) {
    return(strsplit(name, split = "\\.")[[1]][1])
  }))
  
  res <- matrix(rep(0, length(sigs.consensus) * length(tumour.types)),
                ncol = length(tumour.types),
                dimnames = list(sigs.consensus, tumour.types))
  
  bin <- sapply(indv.names, function(name) {
    res[clusters[name], strsplit(name, split = "\\.")[[1]][1]] <<- 1
  })
  
  return(res)
}

# Saves data into an RData file with a customised object name.
# @param  data        data to be saved into an RData file
# @param  file.out    an {@code RDataFileOut} object
# @param  obj.name    the customised variable name of an RData object
Save2RData <- function(data, file.out, obj.name = NULL) {
  UseMethod("Save2RData", file.out)
}

Save2RData.RDataFileOut <- function(data, file.out, obj.name = NULL) {
  obj.name <- ifelse(is.character(obj.name) && trimws(obj.name) != "",
                     trimws(obj.name), file.out$name)
  assign(obj.name, data)
  
  cat("Info: Saving file \"", basename(file.out$fullname), "\" ...\n", sep = "")
  save(list = obj.name, file = file.out$fullname)
  cat("Info: File \"", basename(file.out$fullname), "\" saved\n", sep = "")
}
