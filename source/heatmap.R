library(ggplot2)

if (!exists("Transform3", mode = "function")) {
  source("transform.R")
}

# Public function to plot heatmap for 96 mutation types.
# This function calls the function "Transform3", so the "mut.ctx" shares
# the same constraints of that one. Currenly it accepts only one type of
# "plotter", namely, "HeatmapPlotter"
PlotHeatmap <- function(mut.ctx, plotter, grf.name = "") {
  if (!is.object(plotter) || tail(class(plotter), n = 1) != "HeatmapPlotter") {
    stop("Invalid 'plotter'", call. = F)
  }
  
  cat("Info: Preparing heatmap data ...\n")
  data <- Transform3(mut.ctx, "base.ctx.heatmap.plot")
  grf.name <- ifelse(is.character(grf.name) && trimws(grf.name) != "",
                     trimws(grf.name),
                     strsplit(plotter$file.out$name, "\\.")[[1]][1])
  cat("Info: Heatmap data ready\n")
  
  cat("Info: Start plotting heatmap ...\n")
  ggplot(data, aes(x = nxt_base, y = fwd_base, fill = value)) +
    facet_grid(sample_name ~ mut_base, labeller = label_context) +
    geom_tile() +
    scale_y_discrete(limits = rev(kNtBase)) +
    scale_fill_gradient2(name = "Percentage (log10)", limits= c(-2, 2),
                         low = "white", mid = "gold", high = "red",
                         midpoint = 0,
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = 12)) +
    coord_equal() +
    labs(title = paste("Mutation Heatmap for", grf.name),
         x = "Three Prime Base", y = "Five Prime Base") + 
    plotter$theme
  
  file.fullname <- GenFullFileName(plotter$file.out)
  ggsave(file.fullname, height = max(6, 2.6 + 0.85 * ncol(mut.ctx)),
         units = "in", limitsize = F)
  cat("Info: Plot completed in \"", basename(file.fullname), "\"\n", sep = "")
}
