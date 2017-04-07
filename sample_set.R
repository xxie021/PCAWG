library(reshape2)

if (!exists("VariantSample", mode = "function")) {
  source("variant_sample.R")
}

# Inheritance
# Top: SampleSet
# 2nd Level: SingleSsmSet
SampleSet <- function(sample.set, tumour.name, order.base = "") {
  me <- list(
    samples = sample.set,
    tumour = tumour.name,
    order.base = order.base
  )
  
  class(me) <- append(class(me), "SampleSet")
  return(me)
}

SingleSsmSet <- function(sample.set, tumour.name, order.base = "") {
  class.name <- "SingleSsmSet"
  
  cat("[", class.name, "] Verifying input samples ...\n", sep = "")
  set.validity.marker <- sapply(sample.set, function(ssm) {
    tail(class(ssm), n = 1) == "SingleSSM"
  })
  
  if (anyNA(set.validity.marker) || !all(set.validity.marker, na.rm = T)) {
    stop("Invalid data found in the set", call. = F)
  }
  cat("[", class.name, "] Input samples verified\n", sep = "")
  
  cat("[", class.name, "] Constructing class ...\n", sep = "")
  me <- SampleSet(sample.set, tumour.name, order.base)
  
  cat("[", class.name, "] Creating spectrum data ...\n", sep = "")
  me$data.spectrum <- as.data.frame(sapply(sample.set,
                                           function(ssm) {
                                             ssm$data.spectrum$percentage
                                           }), row.names = kMutBase)
  colnames(me$data.spectrum) <- lapply(sample.set, function(ssm) {
    strsplit(GenerateSampleName(ssm), "-")[[1]][1]
  })
  cat("[", class.name, "] Stacked spectrum data created\n", sep = "")
  
  cat("[", class.name, "] Creating heatmap data ...\n", sep = "")
  list.data.hmap <- lapply(sample.set, function(ssm) ssm$data.hmap)
  me$data.hmap <- do.call("rbind", list.data.hmap)
  cat("[", class.name, "] Heatmap data created\n", sep = "")
  
  class(me) <- append(class(me), class.name)
  cat("[", class.name, "] Class constructed\n", sep = "")
  
  return(me)
}

Summary.SingleSsmSet <- function(obj, type, out.path = NA) {
  type <- tolower(trimws(type))
  
  if (type == "spectrum") {
    smry.data <- t(round(obj$data.spectrum * 100, digits = 2))
    colnames(smry.data) <- c("C:G > A:T", "C:G > G:C", "C:G > T:A",
                             "T:A > A:T", "T:A > C:G", "T:A > G:C")
  } else if (type == "heatmap") {
    smry.data <- round(obj$data.hmap[1:24] * 100, digits = 2)
    row.names(smry.data) <- paste(obj$data.hmap$sample.name, "|",
                                  obj$data.hmap$`5base`)
  } else if (startsWith(type, "heatmap.comp")) {
    smry.data <- do.call("rbind", lapply(obj$samples, function(s) {
      df.ctx <- t(as.data.frame(c(s$c2a.ctx, s$c2g.ctx, s$c2t.ctx,
                                  s$t2a.ctx, s$t2c.ctx, s$t2g.ctx)))
      df.ctx <- cbind(df.ctx, sum(s$c2a.ctx), sum(s$c2g.ctx), sum(s$c2t.ctx),
                      sum(s$t2a.ctx), sum(s$t2c.ctx), sum(s$t2g.ctx))
      if (!endsWith(type, ".num")) {
        df.ctx <- df.ctx / s$total.mut
      }
      row.names(df.ctx) <- GenerateSampleName(s, "short")
      colnames(df.ctx) <- c(sapply(colnames(df.ctx)[1:96], function(x) {
        colname.parts <- unlist(strsplit(x, " |\\."))
        return(paste0(colname.parts[1], "_", colname.parts[2],
                      substr(colname.parts[1], 1, 1), colname.parts[3]))
      }), kMutBase)
      return(df.ctx)
    }))
  } else {
    stop("No data found", call. = F)
  }
  
  cat("Tumour Name:\t", obj$tumour, "\n")
  cat("# Samples:\t", length(obj$samples), "\n")
  print(smry.data)
  cat("\n")
  
  if (!is.na(out.path)) {
    write.table(smry.data,
                file = paste0(out.path, tolower(obj$tumour), ".",
                              type, ".summary.tsv"),
                quote = F, sep = "\t", col.names = NA)
  }
}

PlotSpectrum.SingleSsmSet <- function(obj, plotter, order.desc = F) {
  class.name <- tail(class(obj), n = 1)
  
  cat("[", class.name, "] Start plotting spectrum ...\n", sep = "")
  plot.data <- obj$data.spectrum
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

PlotHeatmap.SingleSsmSet <- function(obj, plotter) {
  class.name <- tail(class(obj), n = 1)
  
  cat("[", class.name, "] Start plotting heatmap ...\n", sep = "")
  plot.data <- melt(obj$data.hmap, id.vars = c("5base", "sample.name"),
                    variable.name = "3base", value.name = "percentage")
  mutbase.plus.3base <- colsplit(plot.data$`3base`, pattern = " ",
                                 names = c("mut.base", "3base"))
  plot.data$mut.base <- mutbase.plus.3base$mut.base
  plot.data$`3base` <- mutbase.plus.3base$`3base`
  plot.data$percentage <- log10(ifelse(plot.data$percentage < 0.0001, 0.01,
                                       plot.data$percentage * 100))
  
  ggplot(plot.data, aes(x = `3base`, y = `5base`, fill = percentage)) +
    facet_grid(sample.name ~ mut.base, labeller = label_context) +
    geom_tile() +
    scale_y_discrete(limits = rev(kNtBase)) +
    scale_fill_gradient2(name = "Percentage (log10)", limits = c(-2, 2),
                         low = "white", mid = "gold", high = "red",
                         midpoint = 0,
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = 12)) +
    coord_equal() +
    labs(title = paste0("Mutation Heatmap for\nCancer \"", obj$tumour, "\""),
         x = "Three Prime Base", y = "Five Prime Base") +
    plotter$theme
  
  ggsave(paste0(plotter$out.path, tolower(obj$tumour),
                ".heatmap.", plotter$type),
         height = max(6, 2.6 + 0.85 * length(obj$samples)), units = "in",
         limitsize = F)
  cat("[", class.name, "] Plotting completed in ", toupper(plotter$type), "\n",
      sep = "")
}

CheckData.SingleSsmSet <- function(obj) {
  print(obj$data.hmap)
}