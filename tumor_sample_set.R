library(reshape2)

if (!exists("VariantSample", mode = "function")) {
  source("variant_sample.R")
}

# Inheritance
# Top: TumorSampleSet
# 2nd Level: TumorSsmSetTypeI
TumorSampleSet <- function(sample.set, tumor.name, order.base = "") {
  me <- list(
    samples = sample.set,
    tumor = tumor.name,
    order.base = order.base
  )
  
  class(me) <- append(class(me), "TumorSampleSet")
  return(me)
}

TumorSsmSetTypeI <- function(sample.set, tumor.name, order.base = "") {
  set.validity.marker <- sapply(sample.set, function(ssm) {
    tail(class(ssm), n = 1) == "SsmTypeI"
  })
  
  if (anyNA(set.validity.marker) || !all(set.validity.marker, na.rm = T)) {
    stop("Invalid data found in the set", call. = F)
  }
  
  me <- TumorSampleSet(sample.set, tumor.name, order.base)
  
  out.data <- as.data.frame(sapply(sample.set,
                                   function(ssm) {
                                     ssm$out.data$percentage
                                   }), row.names = kCategoryTypeI)
  colnames(out.data) <- lapply(sample.set, function(ssm) {
    strsplit(GenFileNamePrefix(ssm), "-")[[1]][1]
  })
  me[["out.data"]] <- out.data
  
  class(me) <- append(class(me), "TumorSsmSetTypeI")
  return(me)
}

Summary.TumorSsmSetTypeI <- function(obj, smry.out.path = NA) {
  smry.data <- t(round(obj$out.data * 100, digits = 2))
  colnames(smry.data) <- c("C:G > A:T", "C:G > G:C", "C:G > T:A",
                           "T:A > A:T", "T:A > C:G", "T:A > G:C")
  
  cat("Tumor Name:\t", obj$tumor, "\n")
  cat("# Samples:\t", length(obj$samples), "\n")
  print(smry.data)
  cat("\n")
  
  if (!is.na(smry.out.path)) {
    write.table(smry.data, file = paste0(smry.out.path,
                                         tolower(obj$tumor), ".summary.tsv"),
                quote = F, sep = "\t", col.names = NA)
  }
}

PlotSignature.TumorSsmSetTypeI <- function(obj, plotter, order.desc = F) {
  cat("Start plotting \"", obj$tumor, "\" SSM Type I signature ...\n", sep = "")
  
  plot.data <- obj$out.data
  if (ncol(plot.data) > 1 && obj$order.base %in% row.names(plot.data)) {
    plot.data <- plot.data[, order(plot.data[obj$order.base, ],
                                   decreasing = order.desc)]
  }
  
  plot.data$mutation.type <- row.names(plot.data)
  plot.data.reshaped <- melt(plot.data, id.vars = c("mutation.type"),
                             variable.name = "sample.name",
                             value.name = "percentage")
  plot.data.reshaped$sample.name <- as.character(plot.data.reshaped$sample.name)
  plot.data.reshaped$sample.name <- factor(plot.data.reshaped$sample.name,
                                           levels = unique(
                                             plot.data.reshaped$sample.name
                                           ))
  
  ggplot(plot.data.reshaped, aes(x = sample.name, y = percentage,
                                 fill = mutation.type)) +
    geom_bar(stat = "identity", colour = "white", width = 1) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::percent) +
    scale_fill_manual(name = "Mutation\nTypes\n", values = kColourTypeI) +
    labs(title = paste0("Spectrum of SSM Type I Signature for\n\"",
                        obj$tumor, "\""), x = "Samples") +
    plotter$plot.theme + plotter$axis.theme + plotter$legend.theme
  
  ggsave(paste0(plotter$out.path, tolower(obj$tumor),
                ".sig.ssm.typeI.", plotter$type),
         width = max(7, round(7 / 50 * length(obj$samples))), units = "in")
  
  cat("Plotting completed in", toupper(plotter$type) ,"\n")
}
