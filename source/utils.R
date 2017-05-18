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
  
  if (any(sapply(mut.ctx.set, nrow) != nrow(mut.ctx.set[[1]]))) {
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

# Merges multiple SSM motif matrices from the given variant sample set.
# @param  sample.set  a list of the {@code SsmSample} objects
# @return             the merged motif matrix
MergeMotifMatrixFromSample <- function(sample.set) {
  if (!is.list(sample.set) || length(sample.set) < 2) {
    stop("Invalid or too few VariantSample set", call. = F)
  }
  
  if (any(sapply(sample.set, function(spl) class(spl)[3] != "SsmSample"))) {
    stop("Invalid VariantSample found in the set", call. = F)
  }
  
  return(MergeMotifMatrix(lapply(sample.set, function(sample) sample$mut.ctx)))
}

# Saves data into an RData file with a customised object name.
# @param  data      data to be saved into an RData file
# @param  file.out  an {@code RDataFileOut} object
# @param  obj.name  the customised variable name of an RData object
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
