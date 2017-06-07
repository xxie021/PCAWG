Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
suppressPackageStartupMessages(library(SomaticSignatures))
suppressPackageStartupMessages(library(NMF))
suppressPackageStartupMessages(library(bioDist))
suppressPackageStartupMessages(library(deconstructSigs))

# Constants
kNumOfGenomes <- c(10, 20, 30, 50, 100, 200)
kMaxNumOfSignatures <- c(3, 5, 7, 10, 15, 20)

# Removes the maximum set of mutation types in a motif matrix that
# togther account for less than or equal to the threshold argument.
# @param  mut.ctx     a motif matrix
# @param  threshold   the cutoff cumulative fraction to filter out minor types
# @param  normalise   the flag indicating if the removal is performed based on
#                     data distribution rather than the mutation counts
# @return             the motif matrix without minor types
RemoveMinorMutationTypes <- function(mut.ctx, threshold = 0.01,
                                     normalise = F) {
  mut.ctx.in <- mut.ctx
  if (normalise) {
    mut.ctx <- sweep(mut.ctx, 2, colSums(mut.ctx), `/`)
  }
  sorted.mut.ctx.pct <- sort(rowSums(mut.ctx) / sum(mut.ctx))
  col.minor <- names(which(cumsum(sorted.mut.ctx.pct) <= threshold))
  reduced.mut.ctx <- as.matrix(mut.ctx.in[
    !(row.names(mut.ctx.in) %in% col.minor),
  ])
  cat("Info: #removed: ", nrow(mut.ctx.in) - nrow(reduced.mut.ctx),
      "; #remaining: ", nrow(reduced.mut.ctx), "\n", sep = "")
  
  if (ncol(mut.ctx.in) == 1) {
    colnames(reduced.mut.ctx) <- colnames(mut.ctx.in)
  }
  return(reduced.mut.ctx)
}

# Applies Monte Carlo bootstrap resampling to avoid overfitting.
# @param  mut.ctx     a motif matrix
# @param  seed        the seed for resampling
# @return             the resampled motif matrix
McBootstrap <- function(mut.ctx, seed = 42) {
  set.seed(seed)
  boot.mut.ctx <- matrix(apply(mut.ctx, 2, function(genome) {
    rmultinom(1, size = sum(genome), prob = genome)
  }), ncol = ncol(mut.ctx))
  colnames(boot.mut.ctx) <- colnames(mut.ctx)
  row.names(boot.mut.ctx) <- row.names(mut.ctx)
  return(boot.mut.ctx)
}

# Estimates the maximum number of possible signatures based on the given number
# of genomes. When deciphering signatures, a range of possible numbers need to
# be tested. This range normally starts from 2 and to such estimated maximum
# number that follows an exponential-fit. The fitted model is built upon six
# paired points {@code kNumOfGenomes} and {@code kMaxNumOfSignatures}.
# @param  n.genomes     the number of genomes for signature construction
# @param  use.raw       the flag indicating if the six paired points are to be
#                       used to adjust results that are predicted by the model
# @return               the estimated maximum number of possible signatures
EstMaxNumOfSignatures <- function(n.genomes, use.raw = FALSE) {
  # If the raw values are used and the number of genomes passed in is
  # one of them, return directly
  if (use.raw && n.genomes %in% kNumOfGenomes) {
    return(kMaxNumOfSignatures[which(kNumOfGenomes == n.genomes)])
  }
  
  fit.data <- data.frame(n_genome = kNumOfGenomes, n_sig = kMaxNumOfSignatures)
  model <- nls(n_sig ~ a + b * log(n_genome), data = fit.data,
               start = list(a = 0, b = 0))
  estimate <- ceiling(predict(model, data.frame(n_genome = n.genomes)))
  
  # If the raw values are used, adjust the estimated number based on adjacent 
  # values. E.g. if the estimated number of 110 genomes is outside [15, 20],
  # it should be adjusted into this range as 110 falls between 100 and 200
  if (use.raw && n.genomes < max(fit.data$n_genome)) {
    if (n.genomes < min(fit.data$n_genome)) {
      estimate <- min(estimate, min(fit.data$n_sig))
    } else {
      lower <- fit.data$n_sig[max(which(fit.data$n_genome < n.genomes))]
      upper <- fit.data$n_sig[min(which(fit.data$n_genome > n.genomes))]
      estimate <- max(lower, min(estimate, upper))
    }
  }
  
  if (estimate <= 0) {
    stop("Estimation can't be done with the given data", call. = F)
  }
  
  return(estimate)
}

# Constructs signatures based on NMF. RSS and CCC are used jointly to determine
# the optimal number of signatures, a.k.a. "k". The cosine similarity is also
# used to evaluate possible results. Two thresholds can be set, namely,
# {@code cosine.low} and {@code cosine.high}. If the cosine similarity between
# any two single signatures is less than {@code cosine.low}, those signatures
# are said to be completedly different. In contrast, if the cosine similarity
# between any two signatures hits {@code cosine.high}, they are considered as
# redundant and hence one of them should be removed. However, if the cosine
# similarity falls into the range between {@code cosine.low} (inclusive) and
# {@code cosine.high} (exclusive), visual inspection is required to detemine
# if the two signatures are identical not. In this case, the function also
# prints a warning message to alert users.
# @param  nmf           an NMF.rank object
# @param  mut.ctx       a motif matrix
# @param  rss.threshold the threshold of percentile for the "k" selection
#                       based on RSS
# @param  seed          the seed to run NMF
# @param  cosine.low    the lower bound of cosine value to indicate different
#                       signatures
# @param  cosine.high   the upper bound of cosine value to indicate redundant
#                       signatures
# @return               the MutationalSignatures object storing deciphered
#                       signatures
ConstSignatures <- function(nmf, mut.ctx,
                            rss.threshold = 0.8, seed = 42,
                            cosine.low = 0.55, cosine.high = 0.85) {
  UseMethod("ConstSignatures", nmf)
}

ConstSignatures.NMF.rank <- function(nmf, mut.ctx,
                                     rss.threshold = 0.8, seed = 42,
                                     cosine.low = 0.55, cosine.high = 0.85) {
  rss.best <- .EstFromRSS(nmf$measures, rss.threshold)
  cat("Info: Best 'k' from RSS: ", rss.best, "\n", sep = "")
  ccc.best <- .EstFromCCC(nmf$measures)
  cat("Info: Best 'k' from CCC: ", ccc.best, "\n", sep = "")
  
  if (rss.best == ccc.best) {
    cat("Info: Estimated number of signatures: ", rss.best, "\n", sep = "")
    sigs <- identifySignatures(mut.ctx, rss.best, nmfDecomposition, seed = seed)
  } else {
    # If the best numbers estimated from RSS and CCC are not equivalent,
    # check the cosine similarity of construsted signatures for each "k"
    # between the boundaries of {@code rss.best} and {@code ccc.best} from
    # large to small values. I.e. if {@code rss.best} returns 5 and
    # {@code ccc.best} returns 2 and 4, check each "k" from 5 to 2
    sel.range <- seq(from = max(rss.best, ccc.best),
                     to = min(rss.best, ccc.best), by = -1)
    for(k in sel.range) {
      sigs <- identifySignatures(mut.ctx, k, nmfDecomposition, seed = seed)
      cosine <- lsa::cosine(signatures(sigs))
      diag(cosine) <- NA
      
      if (all(cosine < cosine.high, na.rm = T)) {
        cat("Info: Estimated number of signatures: ", k, "\n", sep = "")
        if (any(cosine >= cosine.low, na.rm = T)) {
          cat("Warn: When k = ", k, ", the cosine similarities of some",
              " construsted signatures fall into the predefined ambiguous",
              " range. Visual inspection on 'k's from this value to ",
              min(sel.range), " may be required\n", sep = "")
        }
        break
      } else if (k == min(sel.range)) {
        cat("Info: Estimated number of signatures: ", k, "\n",
            "Warn: All 'k's in the range have been tested but some of them",
            " are still highly similar. Visual inspection may be required\n",
            sep = "")
      }
    }
  }
  
  return(sigs)
}

# Adds removed minor mutation types back to the signature matrix with zeros.
# @param  mt.sigs     a signature matrix
# @return             the signature matrix with all 96 mutation types
AddSsm96MinorMutationTypes <- function(mt.sigs) {
  mut.types <- c("CA", "CG", "CT", "TA", "TC", "TG")
  nt.bases <- c("A", "C", "G", "T")
  
  all.sigs <- matrix(rep(0, 96 * ncol(mt.sigs)), ncol = ncol(mt.sigs))
  colnames(all.sigs) <- colnames(mt.sigs)
  row.names(all.sigs) <- sapply(mut.types, function(base) {
    paste0(rep(base, times = 16), " ",
           rep(nt.bases, each = 4), ".", rep(nt.bases, times = 4))
  })
  
  if (!all(row.names(mt.sigs) %in% row.names(all.sigs))) {
    stop("Unrecognised mutation types are detected in signatures", call. = F)
  }
  
  if (nrow(mt.sigs) != 96) {
    all.sigs[row.names(mt.sigs), ] <- mt.sigs
  } else {
    all.sigs <- mt.sigs
  }
  
  return(all.sigs)
}

# Merges multiple signature sets/matrices generated from different tumour types.
# Those signature matrices may have different number of rows due to the removal
# of minor mutation types. In order to process properly, removed minor types
# are added back prior to merging so that each individual signature matrix has
# the equal number of rows (96 rows).
# @param  list.sigs   a named list of the {@code MutationalSignatures} objects
# @return             a merged signature matrix with 96 rows representing
#                     mutation types
MergeSsm96SignaturesFromTumourTypes <- function(list.sigs) {
  if (!is.list(list.sigs)) {
    stop("Invalid MutationalSignatures list", call. = F)
  }
  
  sig.names <- names(list.sigs)
  if (is.null(sig.names) || any(trimws(sig.names) == "")) {
    stop("MutationalSignatures list must have a name for each element",
         call. = F)
  }
  
  sigs <- do.call(cbind, lapply(names(list.sigs), function(name, list.sigs) {
    if (class(list.sigs[[name]])[1] != "MutationalSignatures") {
      stop("Invalid MutationalSignatures is found in the list", call. = F)
    }
    
    mt.sigs <- signatures(list.sigs[[name]])
    colnames(mt.sigs) <- paste0(name, ".", colnames(mt.sigs))
    return(AddSsm96MinorMutationTypes(mt.sigs))
  }, list.sigs = list.sigs))
  
  return(sigs)
}

# Performs agglomerative hierarchical cluster analysis on mutational signatures.
# @param  mt.sigs     a signature matrix
# @param  dist        one of the seven proximity measures:
#                     (1) euclidean
#                     (2) maximum
#                     (3) manhattan
#                     (4) canberra
#                     (5) binary
#                     (6) cosine
#                     (7) correlation
#                     (8) spearman
# @param  method      one (or an unambiguous abbreviation) of the eight
#                     agglomeration methods to be used:
#                     (1) ward.D
#                     (2) ward.D2
#                     (3) single
#                     (4) complete
#                     (5) average (= UPGMA)
#                     (6) mcquitty (= WPGMA)
#                     (7) median (= WPGMC)
#                     (8) centroid (= UPGMC)
# @param  seed        the seed for hierarchical cluster analysis
# @return             the fit of model found by the cluster analysis
ClusterSignatures <- function(mt.sigs, dist = "cosine",
                              method = "complete", seed = 42) {
  if (tolower(dist) == "spearman") {
    distances <- spearman.dist(t(mt.sigs))
  } else {
    distances <- proxy::dist(t(mt.sigs), method = dist)
  }
  set.seed(seed)
  fit <- hclust(distances, method = method)
  
  return(fit)
}

# Summaries a consensus mutational signature from a cluster baesd on
# weighted average of mutation counts of each signature.
# @param  mt.sigs     a signature matrix of a cluster
# @param  sig.name    the string name of the consensus mutational signature
# @return             the summarised consensus mutational signature matrix
GenConsensusSignature <- function(mt.sigs, sig.name = NULL) {
  if (ncol(mt.sigs) == 1) {
    sig.consensus <- mt.sigs
  } else {
    sig.consensus <- as.matrix(
      rowSums(sweep(mt.sigs, 2, colSums(mt.sigs), `*`) / sum(mt.sigs)))
  }
  
  if (is.character(sig.name) && length(sig.name) == 1 &&
      trimws(sig.name) != "") {
    colnames(sig.consensus) <- sig.name
  }
  return(sig.consensus)
}

# Summaries consensus mutational signatures for multiple clusters.
# @param  mt.sigs     an aggregated signature matrix of multiple clusters
# @param  cluster     a named list indicating the cluster indexes for each
#                     tumour type (names)
# @param  sig.names   the list of string names to be assigned to summarised
#                     consensus mutational signatures, used as the prefix if
#                     containing only one name or its length must be equal
#                     to the number of clusters
# @return             the summarised consensus mutational signatures matrix
#                     for each cluster
GenConsensusSignatures <- function(mt.sigs, clusters, sig.names = NULL) {
  if (!identical(sort(colnames(mt.sigs)), sort(names(clusters)))) {
    stop("'mt.sigs' and 'clusters' are NOT related", call. = F)
  }
  
  if (is.null(sig.names) || !is.character(sig.names)) {
    sig.names <- "S"
  }
  
  if (length(sig.names) == 1) {
    sig.names <- paste0(sig.names, ".", sort(unique(clusters)))
  } else if (length(sig.names) != length(unique(clusters))) {
    stop("Length of 'clusters' and 'sig.names' are NOT equal", call. = F)
  }
  
  sigs.consensus <- do.call("cbind", lapply(sort(unique(clusters)), function(idx) {
    sigs.names.single.clust <- names(clusters)[which(clusters == idx)]
    sigs.single.clust <- as.matrix(mt.sigs[, sigs.names.single.clust])
    return(GenConsensusSignature(sigs.single.clust))
  }))
  
  colnames(sigs.consensus) <- sig.names
  return(sigs.consensus)
}

# Evaluates reconstruction error of each sample against multiple consensus
# mutational signatures.
# @param  mt.samples  a sample (96 rows/mutation types) matrix
# @param  mt.sigs     a matrix of consensus mutational signatures, using
#                     COSMIC published signatures if NULL
# @return             the evaluation matrix containing number of included
#                     signatures, reconstruction error and signature weights
#                     for each sample
EvalSsm96Signatures <- function(mt.samples, mt.sigs = NULL) {
  if (is.matrix(mt.samples) && nrow(mt.samples) == 96) {
    df.samples <- as.data.frame(t(mt.samples))
    colnames(df.samples) <- colnames(signatures.cosmic)
  } else {
    stop("Invalid 'mt.samples'. Must be a matrix with 96 rows", call. = F)
  }
  
  if (is.null(mt.sigs)) {
    df.sigs <- signatures.cosmic
  } else if (is.matrix(mt.sigs) && nrow(mt.sigs) == 96) {
    df.sigs <- as.data.frame(t(mt.sigs))
    colnames(df.sigs) <- colnames(signatures.cosmic)
  } else {
    stop("Invalid 'mt.sigs'. Must be a matrix with 96 rows", call. = F)
  }
  
  res <- do.call(rbind, lapply(row.names(df.samples), function(sample.name) {
    cat("Info: Evaluating sample '", sample.name, "' ... \n", sep = "")
    sample <- df.samples[sample.name, ]
    n.muts <- sum(sample)
    stats <- whichSignatures(tumor.ref = sample, signatures.ref = df.sigs,
                             contexts.needed = T)
    n.sigs <- length(which(stats$weights > 0))
    cat("Info: #signatures contributed: ", n.sigs, "\n", sep = "")
    err <- round(sqrt(sum(stats$diff * stats$diff)), digits = 3)
    cat("Info: Error: ", err, "\n", sep = "")
    
    return(c(n.muts, n.sigs, err, stats$weights[1, ]))
  }))
  
  res <- matrix(unlist(res), nrow = nrow(df.samples),
                dimnames = list(
                  row.names(df.samples),
                  c("n.mutations", "n.sigs", "error", row.names(df.sigs))
                ))
  return(res)
}

###################
# Private functions
###################

.EstFromRSS <- function(nmf.meas, threshold = 0.8) {
  if (nrow(nmf.meas) < 3) {
    return(nmf.meas$rank[nrow(nmf.meas)])
  }
  
  diff <- abs(diff(nmf.meas$rss))
  cutoff <- quantile(diff, probs = threshold)
  estimate <- nmf.meas$rank[max(which(diff >= cutoff)) + 1]
  return(estimate)
}

.EstFromCCC <- function(nmf.meas) {
  if (nrow(nmf.meas) == 1) {
    return(nmf.meas$rank[1])
  }
  
  diff <- diff(nmf.meas$cophenetic)
  n.desc <- length(which(diff < 0))
  if (n.desc == 0) {
    # If CCC never decreases, pick up the "k" having the highest value
    estimate <- nmf.meas$rank[nrow(nmf.meas)]
  } else if (n.desc == 1) {
    # If CCC has only one drop, pick up the "k" starting that drop
    estimate <- nmf.meas$rank[which(diff < 0)[1]]
  } else {
    # If CCC has more than one drop, select two "k"s starting the top-two most
    # significant drops and then pick up the one with higher silhouette value
    sorted.diff <- sort(diff)
    estimates <- nmf.meas[which(diff %in% sorted.diff[1:2]), ]
    estimate <- estimates[
      with(estimates, order(silhouette.consensus, decreasing = T)),
    ]$rank[1]
  }
  return(estimate)
}
