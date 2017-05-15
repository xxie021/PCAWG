Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
suppressPackageStartupMessages(library(SomaticSignatures))
suppressPackageStartupMessages(library(NMF))

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
          warning(paste0("When k = ", k, ", the cosine similarities of",
                         " some construsted signatures fall into the",
                         " predefined ambiguous range. Visual inspection on",
                         " 'k's from this value to ", min(sel.range),
                         " may be required"), call. = F)
        }
        break
      } else if (k == min(sel.range)) {
        cat("Info: Estimated number of signatures: ", k, "\n", sep = "")
        warning(paste0("All 'k's in the range have been tested but some of",
                       " them are still highly similar. Visual inspection",
                       " may be required"), call. = F)
      }
    }
  }
  
  return(sigs)
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
