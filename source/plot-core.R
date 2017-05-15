source("source/plot-interface.R")

# Constants
kMutBaseColour <- c("deepskyblue", "black", "tomato",
                    "gray", "yellowgreen", "pink")

PlotSsmCounts.matrix <- function(mut.ctx, plotter, geno.type = "",
                                 log10 = TRUE) {
  data <- reshape2::melt(mut.ctx, varnames = c("type", "id"))
  data$id <- as.factor(data$id)
  if (log10) {
    data$value <- log10(data$value)
  }
  
  geno.type <- ifelse(is.character(geno.type) && trimws(geno.type) != "",
                      trimws(geno.type),
                      strsplit(plotter$file.out$name, "\\.")[[1]][1])
  
  cat("Info: Start plotting mutation type counts ...\n")
  plot <- ggplot(data, aes(x = type, y = value, fill = type)) +
    geom_boxplot() +
    labs(title = paste("Mutation Type Counts of", geno.type),
         x = "Mutation Types",
         y = paste0("Counts", ifelse(log10, " (log10)", ""))) +
    plotter$theme
  
  if (nrow(mut.ctx) > 8) {
    myPalette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
    plot <- plot + scale_fill_manual(values = myPalette(nrow(mut.ctx)))
  } else {
    plot <- plot + scale_fill_brewer(palette = "Dark2")
  }
  
  ggsave(plotter$file.out$fullname, plot = plot, 
         width = max(7, 1.6 + 0.1 * nrow(mut.ctx)), height = 7,
         units = "in", limitsize = F)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotSsmSignatures.MutationalSignatures <- function(ssm.sigs, plotter,
                                                   geno.type = "") {
  if (!is.object(plotter) || class(plotter)[3] != "SignaturePlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  geno.type <- ifelse(is.character(geno.type) && trimws(geno.type) != "",
                      trimws(geno.type),
                      strsplit(plotter$file.out$name, "\\.")[[1]][1])
  
  cat("Info: Start plotting signatures ...\n")
  plot <- plotSignatures(ssm.sigs, normalize = T) + 
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = kMutBaseColour) +
    labs(title = paste("NMF-based Signatures for", geno.type),
         x = "Mutation Types", y = "Percentage of Mutations") +
    plotter$theme
  
  ggsave(plotter$file.out$fullname, plot = plot, 
         width = 7, height = 1.6 + 1.25 * ncol(signatures(ssm.sigs)),
         units = "in", limitsize = F)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotSsmSigContribution.MutationalSignatures <- function(ssm.sigs, plotter,
                                                        geno.type = "",
                                                        order = NULL) {
  if (!is.object(plotter) || class(plotter)[3] != "ContributionPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  contribution <- samples(ssm.sigs)
  contribution <- contribution / rowSums(contribution)
  if (is.character(order) && order %in% colnames(contribution)) {
    row.order <- order(contribution[, which(colnames(contribution) == order)],
                       decreasing = T)
    contribution <- contribution[row.order, ]
  }
  data <- reshape2::melt(contribution, varnames = c("sample", "signature"))
  data$sample = factor(data$sample, levels = unique(data$sample))
  data$signature = factor(data$signature)
  
  geno.type <- ifelse(is.character(geno.type) && trimws(geno.type) != "",
                      trimws(geno.type),
                      strsplit(plotter$file.out$name, "\\.")[[1]][1])
  
  cat("Info: Start plotting signature contribution ...\n")
  plot <- ggplot(data, aes(x = sample, y = value, fill = signature)) +
    geom_bar(stat = "identity", colour = "white", size = 0.8) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = paste("Signature Contribution of", geno.type),
         x = "Samples", y = "Contribution") +
    guides(fill = guide_legend(title = "Signature")) +
    plotter$theme
  
  if (ncol(contribution) > 8) {
    myPalette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
    plot <- plot + scale_fill_manual(values = myPalette(ncol(contribution)))
  } else {
    plot <- plot + scale_fill_brewer(palette = "Set2")
  }
  
  ggsave(plotter$file.out$fullname, plot = plot, 
         width = max(7, 2.5 + 0.15 * nrow(contribution)), height = 7,
         units = "in", limitsize = F)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotCosineSimilarity.MutationalSignatures <- function(ssm.sigs, plotter,
                                                      geno.type = "") {
  if (!is.object(plotter) || class(plotter)[3] != "CosinePlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  cosine <- lsa::cosine(signatures(ssm.sigs))
  cosine[lower.tri(cosine)]<- NA
  data <- reshape2::melt(cosine, na.rm = T)
  label <- round(data$value, 2)
  
  geno.type <- ifelse(is.character(geno.type) && trimws(geno.type) != "",
                      trimws(geno.type),
                      strsplit(plotter$file.out$name, "\\.")[[1]][1])
  
  cat("Info: Start plotting cosine similarity ...\n")
  plot <- ggplot(data, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label)) +
    scale_fill_gradient2(name = "Cosine\n", limit = c(0, 1),
                         low = "white", mid = "turquoise", high = "darkcyan",
                         midpoint = 0.5) +
    coord_fixed() +
    labs(title = paste("Signatures Cosine Similarity of", geno.type)) +
    plotter$theme
  
  ggsave(plotter$file.out$fullname, plot = plot,
         width = 8, height = 7, units = "in")
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotMeasures.NMF.rank <- function(nmf, plotter, geno.type = "") {
  if (!is.object(plotter) || class(plotter)[3] != "MeasurePlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  geno.type <- ifelse(is.character(geno.type) && trimws(geno.type) != "",
                      trimws(geno.type),
                      strsplit(plotter$file.out$name, "\\.")[[1]][1])
  
  cat("Info: Start plotting CCC and RSS measures ...\n")
  plot <- plot(nmf, what = c("cophenetic", "rss"),
               xlab = "Factorization Rank",
               main = paste("NMF Measures of", geno.type)) +
    plotter$theme
  
  ggsave(plotter$file.out$fullname, plot = plot)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}