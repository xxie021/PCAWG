source("source/plot-interface.R")

# Constants
kMutBaseColour <- c("deepskyblue", "black", "tomato",
                    "gray", "yellowgreen", "pink")
kSsimBaseColour <- c(kMutBaseColour, "goldenrod")

PlotSsm6Basics.matrix <- function(mut.ctx, plotter, id, geno.name = NULL,
                                  percentage = FALSE) {
  cat("Info: Preparing basic data ...\n")
  if (is.null(id) || !(id %in% colnames(mut.ctx))) {
    stop("Sample/Genome is not found", call. = F)
  } else {
    id <- as.character(id)
  }
  data <- Transform3(as.matrix(mut.ctx[, id]), "base.summary")
  geno.name <- ifelse(is.character(geno.name) && trimws(geno.name) != "",
                      trimws(geno.name), id)
  
  if (percentage) {
    data <- as.data.frame(data[, "percentage"])
    cat("Info: Basic data ready\n")
    
    cat("Info: Start plotting basic data for sample (in Percentage) ...\n")
    plot <- ggplot(data, aes(x = kMutBase, y = data)) +
      geom_bar(stat = "identity", fill = kMutBaseColour, width = 0.4) +
      scale_y_continuous(breaks = seq(0, max(data), by = 0.05),
                         labels = scales::percent) +
      ggtitle(paste("Percentage of Base Substitutions for", geno.name)) +
      plotter$theme
  } else {
    data <- as.data.frame(data[, "count"])
    cat("Info: Basic data ready\n")
    
    cat("Info: Start plotting basic data for sample (in Count) ...\n")
    plot <- ggplot(data, aes(x = kMutBase, y = data)) +
      geom_bar(stat = "identity", fill = rev(kMutBaseColour), width = 0.4) +
      scale_x_discrete(limits = rev(kMutBase)) +
      scale_y_continuous(labels = scales::comma) +
      coord_flip() +
      labs(title = paste("Number of Base Substitutions for", geno.name),
           y = "Number of Mutations") +
      plotter$theme
  }
  
  ggsave(plotter$file.out$fullname, plot = plot,
         width = 7, height = 7, units = "in")
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotSsm6Spectrum.matrix <- function(mut.ctx, plotter, geno.type = NULL,
                                    order = NULL) {
  cat("Info: Preparing spectrum data ...\n")
  data <- Transform3(mut.ctx, "base.spectrum", order = order)
  geno.type <- ifelse(is.character(geno.type) && trimws(geno.type) != "",
                      trimws(geno.type),
                      strsplit(plotter$file.out$name, "\\.")[[1]][1])
  cat("Info: Spectrum data ready\n")
  
  cat("Info: Start plotting spectrum for set ...\n")
  plot <- ggplot(data, aes(x = sample_name, y = value, fill = mut_base)) +
    geom_bar(stat = "identity", colour = "white", width = 1) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::percent) +
    scale_fill_manual(name = "Mutation\nTypes\n", values = kMutBaseColour) +
    labs(title = paste("Spectrum of Base Substitutions for", geno.type),
         x = "Samples") +
    plotter$theme
  
  ggsave(plotter$file.out$fullname, plot = plot,
         width = max(7, 2.5 + 0.15 * ncol(mut.ctx)), height = 7,
         units = "in", limitsize = F)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotSsm96Heatmap.matrix <- function(mut.ctx, plotter, geno.type = NULL) {
  cat("Info: Preparing heatmap data ...\n")
  data <- Transform3(mut.ctx, "base.ctx.heatmap.plot")
  geno.type <- ifelse(is.character(geno.type) && trimws(geno.type) != "",
                      trimws(geno.type),
                      strsplit(plotter$file.out$name, "\\.")[[1]][1])
  cat("Info: Heatmap data ready\n")
  
  cat("Info: Start plotting heatmap ...\n")
  plot <- ggplot(data, aes(x = nxt_base, y = fwd_base, fill = value)) +
    facet_grid(sample_name ~ mut_base, labeller = label_context) +
    geom_tile() +
    scale_y_discrete(limits = rev(kNtBase)) +
    scale_fill_gradient2(name = "Percentage (log10)", limits= c(-2, 2),
                         low = "white", mid = "gold", high = "red",
                         midpoint = 0,
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = 12)) +
    coord_equal() +
    labs(title = paste("Mutation Heatmap for", geno.type),
         x = "Three Prime Base", y = "Five Prime Base") + 
    plotter$theme
  
  ggsave(plotter$file.out$fullname, plot = plot,
         width = 7, height = max(7, 2.6 + 0.85 * ncol(mut.ctx)),
         units = "in", limitsize = F)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotMutationTypeCounts.matrix <- function(mut.ctx, plotter, geno.type = NULL,
                                          log10 = TRUE) {
  data <- reshape2::melt(mut.ctx, varnames = c("type", "id"))
  data$id <- as.factor(data$id)
  if (log10) {
    data$value <- log10(data$value)
    data$value[which(data$value == -Inf)] <- 0
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
                                                   geno.type = NULL,
                                                   free.yaxis = FALSE) {
  geno.type <- ifelse(is.character(geno.type) && trimws(geno.type) != "",
                      trimws(geno.type),
                      strsplit(plotter$file.out$name, "\\.")[[1]][1])
  
  cat("Info: Start plotting signatures ...\n")
  # Add the prefix "SomaticSignatures::" to distinguish from 
  # the "deconstructSigs" package
  plot <- SomaticSignatures::plotSignatures(ssm.sigs, normalize = T) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = kMutBaseColour) +
    labs(title = paste("NMF-based SSM Signatures for", geno.type),
         x = "Mutation Types", y = "Percentage of Mutations") +
    plotter$theme
  
  if (free.yaxis) {
    plot <- plot + facet_grid(signature ~ alteration, scales = "free_y")
  }
  
  ggsave(plotter$file.out$fullname, plot = plot, 
         width = 7, height = 1.6 + 1.25 * ncol(signatures(ssm.sigs)),
         units = "in", limitsize = F)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotSsim105Signatures.MutationalSignatures <- function(ssim.sigs, plotter,
                                                       geno.type = NULL,
                                                       free.yaxis = FALSE) {
  signatures <- AddSsim105MinorMutationTypes(signatures(ssim.sigs))
  signatures <- sweep(signatures, 2, colSums(signatures), `/`)
  if (sum(signatures["delins", ]) == 0) {
    signatures <- signatures[1:104, ]
  }
  
  data <- reshape2::melt(signatures, varnames = c("full.type", "signature"))
  data$category <- sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", data$full.type)
  data$category <- sub("[ins|del].+", "Indel", data$category)
  data$type <- sub("[ACGTN][ACGTN] (.+)", "\\1", data$full.type)
  data$category <- factor(data$category, levels = unique(data$category))
  data$type <- factor(data$type, levels = unique(data$type))
  
  geno.type <- ifelse(is.character(geno.type) && trimws(geno.type) != "",
                      trimws(geno.type),
                      strsplit(plotter$file.out$name, "\\.")[[1]][1])
  
  cat("Info: Start plotting signatures ...\n")
  plot <- ggplot(data, aes(x = type, y = value, fill = category)) +
    facet_grid(signature ~ category,
               scales = ifelse(free.yaxis, "free", "free_x"),
               space = "free_x") +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = kSsimBaseColour) +
    labs(title = paste("NMF-based SSIM Signatures for", geno.type),
         x = "Mutation Types", y = "Percentage of Mutations") +
    theme_bw() +
    plotter$theme +
    theme(
      axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)
    )
  
  ggsave(plotter$file.out$fullname, plot = plot,
         width = 8, height = 1.6 + 1.25 * ncol(signatures(ssim.sigs)),
         units = "in", limitsize = F)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotSigContribution.MutationalSignatures <- function(sigs, plotter,
                                                     geno.type = NULL,
                                                     order = NULL) {
  contribution <- samples(sigs)
  contribution <- contribution / rowSums(contribution)
  if (is.character(order) && order %in% colnames(contribution)) {
    row.order <- order(contribution[, which(colnames(contribution) == order)],
                       decreasing = T)
    contribution <- contribution[row.order, ]
  } else if (!is.null(order)) {
    cat("Warn: Contribution sorting factor does not exist. Ignored\n")
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

PlotSigPrevalence.matrix <- function(prevalence, plotter, n.samples = NULL) {
  PlotSigPrevalence(as.data.frame(prevalence), plotter, n.samples)
}

PlotSigPrevalence.data.frame <- function(prevalence, plotter,
                                         n.samples = NULL) {
  if (!identical(sort(colnames(prevalence)),
                 c("n.tumour.types", "percentage", "sig"))) {
    stop(paste("Invalid 'prevalence'. Must contain 3 columns:",
               "\"sig\", \"n.tumour.types\" and \"percentage\""), call. = F)
  }
  
  prevalence$sig <- factor(prevalence$sig, levels = prevalence$sig)
  prevalence$percentage[which(prevalence$percentage < 0.01)] <- 0.01
  prevalence$percentage <- (log10(prevalence$percentage) + 2) * 10
  
  title <- "Prevalence of Mutational Signatures"
  if (is.numeric(n.samples) && n.samples > 0) {
    title <- paste0(title, " (#samples = ", as.integer(n.samples), ")")
  }
  
  cat("Info: Start plotting signature prevalence graph ...\n")
  plot <- ggplot(prevalence) +
    geom_bar(aes(x = sig, y = n.tumour.types), stat = "identity",
             colour = "black", fill = "#1B9E77", width = 0.6) +
    geom_line(aes(x = sig, y = percentage, group = 1),
              colour = "#D95F02", size = 0.8) +
    geom_point(aes(x = sig, y = percentage, group = 1),
               colour = "#D95F02", shape= 18, size = 3) +
    geom_text(aes(x = sig, y = percentage,
                  label = paste0(round(10 ^ ((percentage / 10) - 2), 2), "%")),
              colour = "grey30", fontface = "bold", size = 2.5, vjust = -1) +
    scale_x_discrete(expand = c(0.02, 0)) + 
    scale_y_continuous(name = "Number of Cancer Types",
                       breaks = seq(0, 40, by = 2), limits = c(0, 40),
                       sec.axis = sec_axis(
                         trans = ~ . * 2.5,
                         name = "Prevalence in Cancer Samples",
                         labels = c("0.0%", "0.1%", "1.0%",
                                    "10.0%", "100.0%"))) +
    labs(title = title, x = "Signature") +
    plotter$theme
  
  ggsave(plotter$file.out$fullname, plot = plot, 
         width = max(10, 1.5 + 0.2 * nrow(prevalence)), height = 7,
         units = "in", limitsize = F)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotCosineSimilarity.MutationalSignatures <- function(ssm.sigs, plotter,
                                                      geno.type = NULL) {
  PlotCosineSimilarity(signatures(ssm.sigs), plotter, geno.type)
}

PlotCosineSimilarity.matrix <- function(ssm.sigs, plotter, geno.type = NULL) {
  cosine <- lsa::cosine(ssm.sigs)
  cat("Info: Highest cosine value: ", max(cosine[which(cosine < 1)]), "\n",
      sep = "")
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
         width = max(8, 2 + 0.4 * ncol(ssm.sigs)),
         height = max(7, 1.3 + 0.4 * ncol(ssm.sigs)),
         units = "in", limitsize = F)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotMeasures.NMF.rank <- function(nmf, plotter, geno.type = NULL) {
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

PlotLineComparison.list <- function(list.stats, plotter, palette = "Set1",
                                    title = NULL, log10 = FALSE) {
  if (!(length(list.stats) %in% 2:4)) {
    stop("Too few/too many factors. Must be between 2 and 4", call. = F)
  }
  
  if (is.null(names(list.stats)) || is.null(names(list.stats[[1]]))) {
    stop("Invalid factors. Factor or sample names are NOT provided", call. = F)
  }
  
  if (!all(sapply(lapply(list.stats, names),
                  identical, names(list.stats[[1]])))) {
    stop("Factors are NOT relevant", call. = F)
  }
  
  cat("Info: Preparing comparison data ...\n")
  data.raw <- do.call(cbind, list.stats)
  colnames(data.raw) <- names(list.stats)
  if (log10) {
    data.raw <- log10(data.raw)
  }
  data <- reshape2::melt(data.raw)
  data$Var1 <- factor(data$Var1)
  cat("Info: Comparison data ready\n")
  
  cat("Info: Start plotting line comparison graph ...\n")
  plot <- ggplot(data, aes(x = Var1, y = value, colour = Var2, group = Var2)) +
    geom_point(size = 1.2) + 
    geom_line(size = 1) +
    scale_colour_brewer(name = "Factors\n", palette = palette) +
    labs(title = ifelse(is.character(title) && trimws(title) != "",
                        title, "Comparison"),
         x = "Samples") +
    plotter$theme
  
  ggsave(plotter$file.out$fullname, plot = plot,
         width = max(7, 2 + 0.12 * nrow(data.raw)), height = 6,
         units = "in", limitsize = F)
  cat("Info: Plot saved in \"", basename(plotter$file.out$fullname), "\"\n",
      sep = "")
}

PlotDendrogram.hclust <- function(fit, file.out, n.clust = NULL,
                                  group.frame = TRUE, topic = NULL,
                                  ref.line = 0) {
  d <- as.dendrogram(fit)
  
  if (is.null(n.clust)) {
    cat("Info: Number of clusters NOT provided. Estimating ...\n")
    n.clust <- find_k(d, krange = 2:(nleaves(d) - 1))$k
    cat("Info: Estimated number of clusters: ", n.clust, "\n", sep = "")
  } else if (!is.numeric(n.clust) || length(n.clust) > 1 || n.clust < 1) {
    stop("Invalid 'n.clust'. Must be a single positive Integer", call. = F)
  } else {
    n.clust <- as.integer(n.clust)
    cat("Info: Input number of clusters: ", n.clust, "\n", sep = "")
  }
  
  heights <- heights_per_k.dendrogram(d)
  h.max <- ceiling(max(heights))
  h.cut <- heights[as.character(n.clust)]
  cat("Info: Cutting height: ", h.cut, "\n", sep = "")
  
  colours <- RColorBrewer::brewer.pal(8, "Dark2")
  if (n.clust > 8) {
    colours <- rep(colours, times = ceiling(n.clust / 8))[1:n.clust]
  }
  d <- color_labels(d, k = n.clust, col = colours)
  labels_cex(d) <- 0.7
  d <- rev(d)
  
  if (!is.null(topic) && trimws(topic) != "") {
    title <- paste0("Cluster Dendrogram of ", topic, "\n(k = ", n.clust, ")")
  } else {
    title <- paste0("Cluster Dendrogram\n(k = ", n.clust, ")")
  }
  
  cat("Info: Start plotting cluster dendrogram ...\n")
  if (inherits(file.out, "JpgFileOut")) {
    jpeg(filename = file.out$fullname,
         width = 8, height = max(6, 3.5 + 0.11 * nleaves(d)),
         units = "in", res = 300)
  } else {
    pdf(file = file.out$fullname,
        width = 8, height = max(6, 3.5 + 0.11 * nleaves(d)))
  }
  tryCatch({
    #par(bg = "grey92", mar = c(3, 8, 4, 3))
    par(mar = c(3, 8, 4, 3))
    plot_horiz.dendrogram(d, xlim = c(h.max, 0), cex.axis = 0.9)
    title(main = title, line = -3, outer = T)
    if (group.frame) {
      rect.dendrogram(d, k = n.clust, border = "#487AA1", horiz = T,
                      lower_rect = -0.25 * h.max, lwd = 1.1)
    }
    abline(v = h.cut, col = "#F38630", lty = "dashed", lwd = 2)
    if (is.numeric(ref.line) && ref.line > 0 && ref.line < 1) {
      abline(v = ref.line, col = "navy", lty = "dashed", lwd = 2)
    }
  }, finally = {
    dev.off()
  })
  cat("Info: Plot saved in \"", basename(file.out$fullname), "\"\n", sep = "")
}
