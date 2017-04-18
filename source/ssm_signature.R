Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
library(SomaticSignatures)

# Constants
kNumOfGenomes <- c(10, 20, 30, 50, 70, 100, 200)
kMaxNumOfSignatures <- c(3, 5, 7, 10, 10, 15, 20)

# Estimate the maximum number of possible signatures based on the given set of
# genomes. When estimating the number of signatures, a range of possible numbers
# need to be tested. This range normally starts from 2 and to a maximum number
# based on the size of genomes. When the genome size rises to more than 200,
# the exponential-fitted estimation may be far than sufficient and hence causes
# waste on computation. The parameter "use.10pct" can then be set to TRUE
# so that the estimation won't go beyond 10% of the genome size
EstMaxNumOfSignatures <- function(n.genomes, use.10pct = FALSE) {
  # If the number of genomes passed in is one of the known data, return directly
  if (n.genomes %in% kNumOfGenomes) {
    return(kMaxNumOfSignatures[which(kNumOfGenomes == n.genomes)])
  }
  
  fit.data <- data.frame(n_genome = kNumOfGenomes, n_sig = kMaxNumOfSignatures)
  model <- nls(n_sig ~ exp(a + b * n_genome), data = fit.data,
               start = list(a = 0, b = 0))
  estimated <- round(predict(model, list(n_genome = n.genomes)))
  
  # Adjust the estimated number based on adjacent values.
  # E.g. if the estimated number of 110 genomes is outside [15, 20],
  # it should be adjusted into this range as 110 falls between 100 and 200
  if (n.genomes < max(fit.data$n_genome)) {
    if (n.genomes < min(fit.data$n_genome)) {
      estimated <- min(estimated, min(fit.data$n_sig))
    } else {
      lower <- fit.data$n_sig[max(which(fit.data$n_genome < n.genomes))]
      upper <- fit.data$n_sig[min(which(fit.data$n_genome > n.genomes))]
      estimated <- max(lower, min(estimated, upper))
    }
    
  } else if (use.10pct) {
    estimated <- min(n.genomes * 0.1, estimated)
  }
  
  if (estimated > n.genomes) {
    stop("Estimation can't be done with the given data", call. = F)
  }
  
  return(estimated)
}

# Estimate the number of possible signatures based on NMF. RSS is measured
# for the selection of such number
EstNumOfSignatures <- function(mut.ctx, range, n.iter = 500,
                               threshold = 0.8, plot = FALSE) {
  res.raw <- assessNumberSignatures(mut.ctx, range, nReplicates = n.iter)
  res <- dcast(res.raw[, colnames(res.raw) != "ExplainedVariance"],
               NumberSignatures ~ Replicate, value.var = "RSS")
  res <- cbind(res, avg = rowMeans(res[, 2:ncol(res)]))
  res <- cbind(res, diff = c(
    NA, abs(res[2:nrow(res), "avg"] - res[1:(nrow(res) - 1), "avg"])))
  res.cutoff <- quantile(res[, "diff"], probs = threshold, na.rm = T)
  estimated <- res[max(which(res[, "diff"] >= res.cutoff)) + 1,
                   "NumberSignatures"] 
  
  if (plot) {
    print(plotNumberSignatures(res.raw))
  }
  
  return(estimated)
}
