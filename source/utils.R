# Remove the maximum set of any mutation types in a motif matrix that
# togther account for less than or equal to the "threshold"
RemoveMinorMutationTypes <- function(mut.ctx, threshold = 0.01) {
  sorted.mut.ctx.pct <- sort(rowSums(mut.ctx) / sum(mut.ctx))
  col.minor <- names(which(cumsum(sorted.mut.ctx.pct) <= threshold))
  reduced.mut.ctx <- as.matrix(mut.ctx[!(row.names(mut.ctx) %in% col.minor), ])
  if (ncol(mut.ctx) == 1) {
    colnames(reduced.mut.ctx) <- colnames(mut.ctx)
  }
  return(reduced.mut.ctx)
}

# Perform Monte Carlo bootstrap resampling to a motif matrix
McBootstrap <- function(mut.ctx, seed = 42) {
  set.seed(seed)
  boot.mut.ctx <- matrix(apply(mut.ctx, 2, function(genome) {
    rmultinom(1, size = sum(genome), prob = genome)
  }), ncol = ncol(mut.ctx))
  colnames(boot.mut.ctx) <- colnames(mut.ctx)
  row.names(boot.mut.ctx) <- row.names(mut.ctx)
  return(boot.mut.ctx)
}

# Merge multiple motif matrices. If they are to be merged by tumour types,
# i.e. each set of motif matrix represents a specific type of tumour, all
# sample-based counts will be added up prior to merging
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
      stop("No. of motif matrices is not equal to No. of tumour types.",
           call. = F)
    }
    merged <- do.call(cbind, lapply(mut.ctx.set, function(mut.ctx) {
      rowSums(mut.ctx)
    }))
    colnames(merged) <- tumour.types
  } else {
    if (length(mut.ctx.set) < 2) {
      stop("Single motif matrix found. Nothing to merge", call. = F)
    }
    merged <- do.call(cbind, mut.ctx.set)
  }
  
  return(merged)
}

# Merge multiple motif matrices from given variant sample set
MergeMotifMatrixFromSample <- function(sample.set) {
  if (!is.list(sample.set) || length(sample.set) < 2) {
    stop("Invalid or too few VariantSample set", call. = F)
  }
  
  if (any(sapply(sample.set, function(sample) class(sample)[3] !=
                 "SsmSample"))) {
    stop("Invalid VariantSample found in the set", call. = F)
  }
  
  return(MergeMotifMatrix(lapply(sample.set, function(sample) sample$mut.ctx)))
}

# Save data into a RData file with a customised object name
Save2RData <- function(data, file.out, obj.name = NA) {
  if (!is.object(file.out) || class(file.out)[3] != "RDataFileOut") {
    stop("Invalid 'file.out'", call. = F)
  }
  
  if (!exists("GenFullFileName", mode = "function")) {
    source("source/models/file_out.R")
  }
  
  obj.name <- ifelse(is.character(obj.name) && trimws(obj.name) != "",
                     trimws(obj.name), file.out$name)
  assign(obj.name, data)
  file.fullname <- GenFullFileName(file.out)
  
  cat("Info: Saving file \"", basename(file.fullname), "\" ...\n", sep = "")
  save(list = obj.name, file = file.fullname)
  cat("Info: File \"", basename(file.fullname), "\" saved\n", sep = "")
}
