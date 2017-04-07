if (!exists("VariantSample", mode = "function")) {
  source("variant_sample.R")
}

# Inheritance
# Top: VariantSet
# 2nd Level: SsmSampleSet, SsmTumourSet
VariantSet <- function(variant.set, set.class, order.base = "") {
  class.name <- "VariantSet"
  
  cat("[", class.name, "] Verifying input variant set ...\n", sep = "")
  if (!VerifyVariantSet(variant.set, set.class)) {
    stop("Invalid data found in the input dataset", call. = F)
  }
  cat("[", class.name, "] Input variant set verified\n", sep = "")
  
  me <- list(
    variant.set = variant.set,
    order.base = order.base
  )
  
  class(me) <- append(class(me), class.name)
  return(me)
}

SsmSampleSet <- function(sample.set, tumour.name, order.base = "") {
  class.name <- "SsmSampleSet"
  
  cat("[", class.name, "] Constructing class ...\n", sep = "")
  me <- VariantSet(sample.set, "SingleSSM", order.base)
  me$tumour <- tumour.name
  
  cat("[", class.name, "] Creating spectrum data on base substitutions ...\n",
      sep = "")
  me$mut.base <- sapply(sample.set, function(x) {
    x$mut.base[, "percentage"]
  })
  colnames(me$mut.base) <- lapply(sample.set, function(x) {
    strsplit(GenerateSampleName(x), "-")[[1]][1]
  })
  cat("[", class.name, "] Spectrum data on base substitutions created\n",
      sep = "")
  
  cat("[", class.name, "] Creating spectrum data with context ...\n", sep = "")
  me$mut.ctx <- do.call("cbind", lapply(sample.set, function(x) x$mut.ctx))
  me$mut.ctx.aggr <- as.matrix(rowSums(me$mut.ctx))
  colnames(me$mut.ctx.aggr) <- tumour.name
  me$mut.total <- sum(me$mut.ctx.aggr)
  cat("[", class.name, "] Spectrum data with context created\n", sep = "")
  
  class(me) <- append(class(me), class.name)
  cat("[", class.name, "] Class constructed\n", sep = "")
  
  return(me)
}

SsmTumourSet <- function(tumour.set, order.base = "") {
  class.name <- "SsmTumourSet"
  
  cat("[", class.name, "] Constructing class ...\n", sep = "")
  me <- VariantSet(tumour.set, "SsmSampleSet", order.base)
  
  cat("[", class.name, "] Creating spectrum data with context ...\n", sep = "")
  me$mut.ctx <- do.call("cbind", lapply(tumour.set, function(x) x$mut.ctx.aggr))
  me$mut.total <- sum(me$mut.ctx)
  cat("[", class.name, "] Spectrum data with context created\n", sep = "")
  
  class(me) <- append(class(me), class.name)
  cat("[", class.name, "] Class constructed\n", sep = "")
  
  return(me)
}

VerifyVariantSet <- function(dataset, class.exp) {
  set.validity.marker <- sapply(dataset, function(x) {
    tail(class(x), n = 1) == class.exp
  })
  
  if (anyNA(set.validity.marker) || !all(set.validity.marker, na.rm = T)) {
    return(FALSE)
  }
  
  return(TRUE)
}

CheckMinorMutationTypes <- function(obj, threshold = "0.01") {
  if (!"VariantSet" %in% class(obj)) {
    stop("'obj' is not a type or a subtype of VariantSet", call. = F)
  }
  
  return(names(which(rowSums(obj$mut.ctx) / obj$mut.total <= threshold)))
}

Summary.VariantSet <- function(obj, type, out.path = NA) {
  type <- tolower(trimws(type))
  set.class <- tail(class(obj), n = 1)
  
  if (type == "mut.base" && set.class == "SsmSampleSet") {
    smry.data <- t(round(obj$mut.base * 100, digits = 2))
    colnames(smry.data) <- kMutBase2
  } else if (startsWith(type, "mut.ctx")) {
    smry.data <- t(obj$mut.ctx)
    if (endsWith(type, ".pct")) {
      smry.data <- sweep(smry.data, 2, rowSums(smry.data), `/`)
    }
    smry.data <- cbind(smry.data, rowSums(smry.data[, 1:16]),
                       rowSums(smry.data[, 17:32]),
                       rowSums(smry.data[, 33:48]),
                       rowSums(smry.data[, 49:64]),
                       rowSums(smry.data[, 65:80]),
                       rowSums(smry.data[, 81:96]))
    colnames(smry.data) <- c(sapply(colnames(smry.data)[1:96], function(x) {
      colname.parts <- unlist(strsplit(x, " |\\."))
      return(paste0(colname.parts[1], "_", colname.parts[2],
                    substr(colname.parts[1], 1, 1), colname.parts[3]))
    }), kMutBase)
  } else if (startsWith(type, "heatmap")) {
    smry.data <- obj$mut.ctx
    if (endsWith(type, ".pct")) {
      smry.data <- round(sweep(smry.data, 2, colSums(smry.data), `/`) * 100,
                         digits = 2)
    }
    sample.names <- colnames(smry.data)
    num.of.samples <- length(sample.names)
    smry.data <- as.data.frame(matrix(smry.data, ncol = 1, byrow = T))
    colnames(smry.data) <- "count"
    smry.data$`5base` <- rep(rep(kNtBase, each = 4), 6 * num.of.samples)
    smry.data$`3base` <- rep(kNtBase, 4 * 6 * num.of.samples)
    smry.data$mut.base <- rep(rep(kMutBase, each = 16), num.of.samples)
    smry.data$sample.name <- rep(sample.names, each = 4 * 4 * 6)
    smry.data <- dcast(smry.data, sample.name + `5base` ~ mut.base + `3base`,
                       value.var = "count")
    row.names(smry.data) <- paste(ttt$sample.name, "|", ttt$`5base`)
    smry.data <- smry.data[, c(-2, -1)]
  } else {
    stop("No data found", call. = F)
  }
  
  cat(rep("#", 20), "\n", sep = "")
  if (set.class == "SsmSampleSet") {
    cat("# Tumour Name:\t", obj$tumour, "\n")
  }
  cat("# No. of Sets:\t", length(obj$variant.set), "\n")
  cat(rep("#", 20), "\n", sep = "")
  print(smry.data)
  cat("\n")
  
  if (!is.na(out.path)) {
    write.table(smry.data,
                file = paste0(out.path, tolower(obj$tumour), ".",
                              type, ".summary.tsv"),
                quote = F, sep = "\t", col.names = NA)
  }
}

PlotBaseSpectrum.SsmSampleSet <- function(obj, plotter, order.desc = F) {
  class.name <- tail(class(obj), n = 1)
  
  cat("[", class.name, "] Start plotting spectrum ...\n", sep = "")
  plot.data <- as.data.frame(obj$mut.base)
  row.names(plot.data) <- kMutBase
  if (ncol(plot.data) > 1 && obj$order.base %in% row.names(plot.data)) {
    plot.data <- plot.data[, order(plot.data[obj$order.base, ],
                                   decreasing = order.desc)]
  }
  
  plot.data$mutation.type <- row.names(plot.data)
  plot.data <- melt(plot.data, id.vars = c("mutation.type"),
                    variable.name = "sample.name", value.name = "percentage")
  plot.data$sample.name <- as.character(plot.data$sample.name)
  plot.data$sample.name <- factor(plot.data$sample.name,
                                  levels = unique(plot.data$sample.name))
  
  ggplot(plot.data, aes(x = sample.name, y = percentage,
                        fill = mutation.type)) +
    geom_bar(stat = "identity", colour = "white", width = 1) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::percent) +
    scale_fill_manual(name = "Mutation\nTypes\n", values = kMutBaseColour) +
    labs(title = paste0("Spectrum of Base Substitutions for\nCancer \"",
                        obj$tumour, "\""), x = "Samples") +
    plotter$theme
  
  ggsave(paste0(plotter$out.path, tolower(obj$tumour),
                ".spectrum.pct.", plotter$type),
         width = max(7, round(7 / 50 * length(obj$samples))), units = "in",
         limitsize = F)
  cat("[", class.name, "] Plotting completed in ", toupper(plotter$type), "\n",
      sep = "")
}

PlotHeatmap.VariantSet <- function(obj, plotter) {
  class.name <- tail(class(obj), n = 1)
  
  cat("[", class.name, "] Preparing heatmap data ...\n", sep = "")
  sample.names <- colnames(obj$mut.ctx)
  num.of.samples <- length(sample.names)
  plot.data <- as.data.frame(matrix(sweep(obj$mut.ctx, 2,
                                          colSums(obj$mut.ctx), `/`),
                                    ncol = 1))
  colnames(plot.data) <- "percentage"
  plot.data$mut.base <- rep(rep(kMutBase, each = 16), times = num.of.samples)
  plot.data$`5base` <- rep(rep(kNtBase, each = 4), times = 6 * num.of.samples)
  plot.data$`3base` <- rep(kNtBase, times = 4 * 6 * num.of.samples)
  plot.data$sample.name <- rep(sample.names, each = 96)
  plot.data$percentage <- log10(ifelse(plot.data$percentage < 0.0001, 0.01,
                                       plot.data$percentage * 100))
  
  if (class.name == "SsmTumourSet") {
    out.name <- "wg"
    plot.title <- "Mutation Heatmap for Whole-Genome Data"
  } else {
    out.name <- obj$tumour
    plot.title <- paste0("Mutation Heatmap for\nCancer \"", obj$tumour, "\"")
  }
  cat("[", class.name, "] Heatmap data ready\n", sep = "")
  
  cat("[", class.name, "] Start plotting heatmap ...\n", sep = "")
  ggplot(plot.data, aes(x = `3base`, y = `5base`, fill = percentage)) +
    facet_grid(sample.name ~ mut.base, labeller = label_context) +
    geom_tile() +
    scale_y_discrete(limits = rev(kNtBase)) +
    scale_fill_gradient2(name = "Percentage (log10)", limits= c(-2, 2),
                         low = "white", mid = "gold", high = "red",
                         midpoint = 0,
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = 12)) +
    coord_equal() +
    labs(title = plot.title, x = "Three Prime Base", y = "Five Prime Base") + 
    plotter$theme
  
  ggsave(paste0(plotter$out.path, tolower(out.name), ".heatmap.", plotter$type),
         height = max(6, 2.6 + 0.85 * length(obj$samples)), units = "in",
         limitsize = F)
  cat("[", class.name, "] Plotting completed in ", toupper(plotter$type), "\n",
      sep = "")
}
