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
MCBootstrap <- function(mut.ctx) {
  set.seed(42)
  boot.mut.ctx <- matrix(apply(mut.ctx, 2, function(genome) {
    rmultinom(1, size = sum(genome), prob = genome)
  }), ncol = ncol(mut.ctx))
  colnames(boot.mut.ctx) <- colnames(mut.ctx)
  row.names(boot.mut.ctx) <- row.names(mut.ctx)
  return(boot.mut.ctx)
}
