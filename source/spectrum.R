library(ggplot2)

if (!exists("Transform3", mode = "function")) {
  source("source/transform.R")
}

# Constants
kMutBaseColour <- c("deepskyblue", "black", "tomato",
                    "gray", "yellowgreen", "pink")

# Public function to plot spectrum for 6 base mutation types.
# This function calls the function "Transform3", so the "mut.ctx" shares
# the same constraints of that one
# Optional parameters:
#   (1) percentage = TRUE/FALSE
#   (2) order (one of "kMutBase")
PlotBaseSpectrum <- function(mut.ctx, plotter, grf.name = "", ...) {
  data.type <- ifelse(ncol(mut.ctx) == 1, "sample", "sample.set")
  args <- list(...)
  
  switch(data.type,
         sample = {
           cat("Info: Single sample detacted\n")
           if ("percentage" %in% names(args)) {
             .PlotBaseSpectrumForSample(mut.ctx, plotter, grf.name,
                                        args$percentage)
           } else {
             .PlotBaseSpectrumForSample(mut.ctx, plotter, grf.name)
           }
         },
         sample.set = {
           cat("Info: Multiple samples detacted\n")
           if ("order" %in% names(args)) {
             .PlotBaseSpectrumForSet(mut.ctx, plotter, grf.name, args$order)
           } else {
             .PlotBaseSpectrumForSet(mut.ctx, plotter, grf.name)
           }
         },
         stop("Unknown error", call. = F)
  )
}

###################
# Private functions
###################

# Plot spectrum of 6 base mutation types for a sample.
# I.e. the motif matrix contains only one column
.PlotBaseSpectrumForSample <- function(mut.ctx, plotter, grf.name = "",
                                       percentage = FALSE) {
  if (!is.object(plotter) || class(plotter)[2] != "XPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  cat("Info: Preparing spectrum data for sample ...\n")
  data <- Transform3(mut.ctx, "base.summary")
  grf.name <- ifelse(is.character(grf.name) && trimws(grf.name) != "",
                     trimws(grf.name),
                     strsplit(plotter$file.out$name, "\\.")[[1]][1])
  
  if (percentage) {
    data <- as.data.frame(data[, "percentage"])
    cat("Info: Spectrum data ready\n")
    
    cat("Info: Start plotting spectrum for sample (in Percentage) ...\n")
    ggplot(data, aes(x = kMutBase, y = data)) +
      geom_bar(stat = "identity", fill = kMutBaseColour, width = 0.4) +
      scale_y_continuous(breaks = seq(0, max(data), by = 0.05),
                         labels = scales::percent) +
      ggtitle(paste0("Percentage of Base Substitutions for \nSample \"",
                     grf.name, "\"")) +
      plotter$theme
  } else {
    data <- as.data.frame(data[, "count"])
    cat("Info: Spectrum data ready\n")
    
    cat("Info: Start plotting spectrum for sample (in Number) ...\n")
    ggplot(data, aes(x = kMutBase, y = data)) +
      geom_bar(stat = "identity", fill = rev(kMutBaseColour), width = 0.4) +
      scale_x_discrete(limits = rev(kMutBase)) +
      scale_y_continuous(labels = scales::comma) +
      coord_flip() +
      labs(title = paste0("Number of Base Substitutions for\nCancer \"",
                          grf.name, "\""), y = "Number of Mutations") +
      plotter$theme
  }
  
  file.fullname <- GenFullFileName(plotter$file.out)
  ggsave(file.fullname)
  cat("Info: Plot completed in \"", basename(file.fullname), "\"\n", sep = "")
}

# Plot spectrum of 6 base mutation types for multiple samples.
# I.e. the motif matrix contains more than one column
.PlotBaseSpectrumForSet <- function(mut.ctx, plbotter, grf.name = "",
                                    order = NA) {
  if (!is.object(plotter) || class(plotter)[3] != "StackedBarPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  cat("Info: Preparing spectrum data for set ...\n")
  data <- Transform3(mut.ctx, "base.spectrum", order = order)
  grf.name <- ifelse(is.character(grf.name) && trimws(grf.name) != "",
                     trimws(grf.name),
                     strsplit(plotter$file.out$name, "\\.")[[1]][1])
  cat("Info: Spectrum data ready\n")
  
  cat("Info: Start plotting spectrum for set ...\n")
  ggplot(data, aes(x = sample_name, y = value, fill = mut_base)) +
    geom_bar(stat = "identity", colour = "white", width = 1) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::percent) +
    scale_fill_manual(name = "Mutation\nTypes\n", values = kMutBaseColour) +
    labs(title = paste0("Spectrum of Base Substitutions for\nCancer \"",
                        grf.name, "\""), x = "Samples") +
    plotter$theme
  
  file.fullname <- GenFullFileName(plotter$file.out)
  ggsave(file.fullname, width = max(7, round(7 / 50 * ncol(mut.ctx))),
         units = "in", limitsize = F)
  cat("Info: Plot completed in \"", basename(file.fullname), "\"\n", sep = "")
}
