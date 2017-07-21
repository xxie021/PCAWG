################################################################################
# NOTE:
# All functions in this script file support context motifs with length = 3 only
# I.e. base subsitutions with their immidate flanking 5' and 3' bases
# DON'T use them for context motifs with other lengths
################################################################################

# Constants
kNtBase <- c("A", "C", "G", "T")
kMutBase <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
kNumOfMutCtx3 <- 96

#' Transforms a motif matrix to the requested format. The ending "3" in the
#' function name indicates the support of context motifs with length = 3 only
#' @param  mut.ctx        a motif matrix (3-mer context, 96 rows)
#' @param  type           one of the five types related to internal functions:
#'                        (1) base.summary
#'                        (2) base.spectrum
#'                        (3) base.ctx.summary
#'                        (4) base.ctx.heatmap.plot
#'                        (5) base.ctx.heatmap.real
#' @param  ...            optional parameters
#'         normalise      the flag indicating if data need to be normalised;
#'                        used with the type "base.summary" only
#'         order          one of the base substitution in \code{kMutBase} to
#'                        order the result; used with the type "base.spectrum"
#'                        only
#'         data.type      either "frac" or "log" indicating if data need to be
#'                        converted into fraction or Log 10 scale; used with
#'                        the type "base.ctx.heatmap.real" only
#' @return                the requested format of the motif matrix
Transform3 <- function(mut.ctx, type, ...) {
  if (nrow(mut.ctx) != kNumOfMutCtx3) {
    stop("Invalid motif matrix. Must have 96 rows", call. = F)
  }
  
  if (!is.character(type)) {
    stop("Invalid data type of 'type'. Must be string", call. = F)
  }
  
  type <- tolower(type)
  args <- list(...)
  
  return(switch(type,
                base.summary = {
                  if ("normalise" %in% names(args)) {
                    .BaseTransform2summary(mut.ctx, args$normalise)
                  } else {
                    .BaseTransform2summary(mut.ctx)
                  }
                },
                base.spectrum = {
                  if ("order" %in% names(args)) {
                    .BaseTransform2spectrum(mut.ctx, args$order)
                  } else {
                    .BaseTransform2spectrum(mut.ctx)
                  }
                },
                base.ctx.summary = .BaseContextTransform2summary(mut.ctx),
                base.ctx.heatmap.plot = 
                  .BaseContextTransform2heatmapPlot(mut.ctx),
                base.ctx.heatmap.real = {
                  if ("data.type" %in% names(args)) {
                    .BaseContextTransform2heatmapReal(mut.ctx, args$data.type)
                  } else {
                    .BaseContextTransform2heatmapReal(mut.ctx)
                  }
                },
                stop("Unrecognised transform type", call. = F)
  ))
}

###################
# Private functions
###################

# Transforms a motif matrix to show info of 6 base mutation types
# of each genome.
.BaseTransform2summary <- function(mut.ctx, normalise = FALSE) {
  data <- cbind(colSums(as.matrix(mut.ctx[1:16, ])),
                colSums(as.matrix(mut.ctx[17:32, ])),
                colSums(as.matrix(mut.ctx[33:48, ])),
                colSums(as.matrix(mut.ctx[49:64, ])),
                colSums(as.matrix(mut.ctx[65:80, ])),
                colSums(as.matrix(mut.ctx[81:96, ])))
  colnames(data) <- kMutBase
  
  if (nrow(data) == 1 && !normalise) {
    data <- t(rbind(data, data / sum(data)))
    colnames(data) <- c("count", "percentage")
  } else if (normalise) {
    data <- data / rowSums(data)
    if (nrow(data) == 1) {
      row.names(data) <- colnames(mut.ctx)
    }
  }
  
  return(data)
}

# Transforms a motif matrix for stacked bar plot on 6 base mutation types
# using ggplot2.
.BaseTransform2spectrum <- function(mut.ctx, order = NULL) {
  data <- t(.BaseTransform2summary(mut.ctx, normalise = T))
  if (ncol(data) > 1) {
    if (is.character(order) && order %in% kMutBase) {
      data <- data[, order(data[order, ], decreasing = T)]
    } else if (!is.null(order)) {
      cat("Warn: Sorting factor does not exist. Ignored\n")
    }
  }
  
  sample.names <- colnames(data)
  n.samples <- length(sample.names)
  
  data <- as.data.frame(matrix(data, ncol = 1))
  colnames(data) <- "value"
  data$mut_base <- rep(kMutBase, times = n.samples)
  data$sample_name <- factor(rep(sample.names, each = length(kMutBase)),
                             levels = sample.names)
  return(data)
}

# Transforms a motif matrix to show counts of 96 mutation types of each genome.
.BaseContextTransform2summary <- function(mut.ctx) {
  data <- t(mut.ctx)
  data <- cbind(data, rowSums(data[, 1:16]), rowSums(data[, 17:32]),
                rowSums(data[, 33:48]), rowSums(data[, 49:64]),
                rowSums(data[, 65:80]), rowSums(data[, 81:96]))
  colnames(data) <- c(sapply(colnames(data)[1:96], function(n) {
    colname.parts <- unlist(strsplit(n, " |\\."))
    return(paste0(colname.parts[1], "_", colname.parts[2],
                  substr(colname.parts[1], 1, 1), colname.parts[3]))
  }), kMutBase)
  return(data)
}

# Transforms a motif matrix for heatmap plot on 96 mutation types using ggplot2.
.BaseContextTransform2heatmapPlot <- function(mut.ctx, log10 = TRUE) {
  data <- sweep(mut.ctx, 2, colSums(mut.ctx), `/`)
  if (log10) {
    data <- log10(ifelse(data < 0.0001, 0.01, data * 100))
  }
  
  data <- .BaseContextTransform2longformat(data)
  return(data)
}

# Transforms a motif matrix to match heatmap representation.
.BaseContextTransform2heatmapReal <- function(mut.ctx, data.type = "frac") {
  data <- mut.ctx
  if (data.type %in% c("frac", "log")) {
    data <- sweep(data, 2, colSums(data), `/`)
  }
  if (data.type == "log") {
    data <- log10(ifelse(data < 0.0001, 0.01, data * 100))
  }
  
  data <- .BaseContextTransform2longformat(data)
  data <- reshape2::dcast(data, sample_name + fwd_base ~ mut_base + nxt_base,
                          value.var = "value")
  row.names(data) <- paste(data$sample_name, "|", data$fwd_base)
  data <- data[, c(-2, -1)]
  return(data)
}

# Transforms a motif matrix to long-format data for 96 mutation types.
.BaseContextTransform2longformat <- function(mut.ctx) {
  sample.names <- colnames(mut.ctx)
  n.samples <- length(sample.names)
  
  data <- as.data.frame(matrix(mut.ctx, ncol = 1))
  colnames(data) <- "value"
  data$mut_base <- rep(rep(kMutBase, each = 4 * 4), times = n.samples)
  data$fwd_base <- rep(rep(kNtBase, each = 4), times = 6 * n.samples)
  data$nxt_base <- rep(kNtBase, times = 6 * 4 * n.samples)
  data$sample_name <- rep(sample.names, each = kNumOfMutCtx3)
  return(data)
}
