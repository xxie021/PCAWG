source("source/plot-interface.R")

# Constants
kMutBaseColour <- c("deepskyblue", "black", "tomato",
                    "gray", "yellowgreen", "pink")

PlotSsm6Basics.matrix <- function(mut.ctx, plotter, id, geno.name = "",
                                  percentage = FALSE) {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
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

PlotSsm6Spectrum.matrix <- function(mut.ctx, plotter, geno.type = "",
                                    order = NULL) {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
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

PlotSsm96Heatmap.matrix <- function(mut.ctx, plotter, geno.type = "") {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
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

PlotSsmCounts.matrix <- function(mut.ctx, plotter, geno.type = "",
                                 log10 = TRUE) {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
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
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
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
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  contribution <- samples(ssm.sigs)
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

PlotCosineSimilarity.MutationalSignatures <- function(ssm.sigs, plotter,
                                                      geno.type = "") {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
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
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
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
