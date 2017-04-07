Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
library(SomaticSignatures)

# Constants
kNtBase <- c("A", "C", "G", "T")
kMutBase <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
kMutBase2 <- c("C:G > A:T", "C:G > G:C", "C:G > T:A",
               "T:A > A:T", "T:A > C:G", "T:A > G:C")

kMutBaseColour <- c("deepskyblue", "black", "tomato",
                    "gray", "yellowgreen", "pink")

# Inheritance
# Top: VariantSample
# 2nd Level: SingleSSM
VariantSample <- function(vcf, ref, id = "", name = "") {
  if (class(vcf) != "CollapsedVCF") {
    stop("Invalid VCF data found", call. = F)
  }
  
  if (class(ref) != "FaFile") {
    stop("Invalid reference sequence found", call. = F)
  }
  
  #tryCatch({
  #  countFa(ref)
  #}, error = function(e) {
  #  stop("Invalid reference sequence found", call. = F)
  #})
  
  if (!is.character(id)) {
    warning("'id' is not type of character. Converted", call. = F)
    id <- as.character(id)
  }
  
  if (!is.character(name)) {
    warning("'name' is not type of character. Converted", call. = F)
    name <- as.character(name)
  }
  
  if (trimws(id) == "" && trimws(name) == "") {
    stop("'id' and 'name' can't be both empty", call. = F)
  }
  
  me <- list(
    vcf = vcf,
    ref = ref,
    id = id,
    name = name
  )
  
  class(me) <- append(class(me), "VariantSample")
  return(me)
}

SingleSSM <- function(vcf, ref, id = "", name = "") {
  class.name <- "SingleSSM"
  
  cat("[", class.name, "] Constructing class ...\n", sep = "")
  me <- VariantSample(vcf, ref, id, name)
  idx <- which(isSNV(me$vcf))
  
  seq.ctx <- mutationContext(as(vcf[idx], "VRanges"), ref)
  sampleNames(seq.ctx) <- id
  motif.num <- motifMatrix(seq.ctx, normalize = F)
  
  me$idx <- idx
  me$mut.total <- sum(motif.num)
  me$c2a.ctx <- CountBaseSubstitutionWithContext(motif.num, "c", "a")
  me$c2g.ctx <- CountBaseSubstitutionWithContext(motif.num, "c", "g")
  me$c2t.ctx <- CountBaseSubstitutionWithContext(motif.num, "c", "t")
  me$t2a.ctx <- CountBaseSubstitutionWithContext(motif.num, "t", "a")
  me$t2c.ctx <- CountBaseSubstitutionWithContext(motif.num, "t", "c")
  me$t2g.ctx <- CountBaseSubstitutionWithContext(motif.num, "t", "g")
  
  cat("[", class.name, "] Creating spectrum data on base substitutions...\n",
      sep = "")
  me$mut.base <- matrix(c(sum(me$c2a.ctx), sum(me$c2g.ctx), sum(me$c2t.ctx),
                          sum(me$t2a.ctx), sum(me$t2c.ctx), sum(me$t2g.ctx)),
                        ncol = 1)
  me$mut.base <- cbind(me$mut.base, me$mut.base / me$mut.total)
  colnames(me$mut.base) <- c("count", "percentage")
  cat("[", class.name, "] Spectrum data on base substitutions created\n",
      sep = "")
  
  cat("[", class.name, "] Creating spectrum data with context...\n", sep = "")
  me$mut.ctx <- motif.num
  cat("[", class.name, "] Spectrum data with context created\n", sep = "")
  
  class(me) <- append(class(me), class.name)
  cat("[", class.name, "] Class constructed\n", sep = "")
  
  return(me)
}

SingleSIM <- function(vcf, ref = "", id = "", name = "") {
  class.name <- "SingleSIM"
  
  cat("[", class.name, "] Constructing class ...\n", sep = "")
  me <- VariantSample(vcf, ref, id, name)
  idx.ins <- which(isInsertion(me$vcf))
  idx.del <- which(isDeletion(me$vcf))
  
  me$idx.ins <- idx.ins
  me$idx.del <- idx.del
  
  class(me) <- append(class(me), class.name)
  cat("[", class.name, "] Class constructed\n", sep = "")
  
  return(me)
}

CountBaseSubstitutionWithContext <- function(obj, nt.ref, nt.alt) {
  if (class(obj) != "matrix") {
    stop("'obj' is not a matrix", call. = F)
  }
  
  mut.base <- toupper(paste0(nt.ref, nt.alt))
  mut.matched.ctx <- obj[which(sapply(strsplit(row.names(obj), " "),
                                      function(x) x[1]) == mut.base), ]
  
  return(mut.matched.ctx)
}

GenerateSampleName <- function(obj, out.type = NA) {
  if (!"VariantSample" %in% class(obj)) {
    stop("'obj' is not a type or a subtype of VariantSample", call. = F)
  }
  
  out.type <- ifelse(is.na(out.type), "normal", tolower(trimws(out.type)))
  
  if (out.type == "full" && obj$id != "" && obj$name != "") {
    return(paste0(obj$id, ".", obj$name))
  }
  
  normal.name <- ifelse(obj$id == "", obj$name, obj$id)
  return(ifelse(out.type == "short", substr(normal.name, 1, 8), normal.name))
}

Summary.SingleSSM <- function(obj, out.path = NA) {
  smry.data <- as.data.frame(obj$mut.base)
  smry.data$percentage <- round(smry.data$percentage * 100, digits = 2)
  row.names(smry.data) <- kMutBase2
  
  cat(rep("#", 20), "\n", sep = "")
  cat("# Sample ID:\t", obj$id, "\n")
  cat("# Sample name:\t", obj$name, "\n")
  cat(rep("#", 20), "\n", sep = "")
  print(smry.data)
  cat("\n")
  
  if (!is.na(out.path)) {
    write.table(smry.data,
                file = paste0(out.path, tolower(GenerateSampleName(obj)),
                              ".summary.tsv"),
                quote = F, sep = "\t", col.names = NA)
  }
}

PlotBaseSpectrum.SingleSSM <- function(obj, plotter, chart.type = "pct") {
  class.name <- tail(class(obj), n = 1)

  cat("[", class.name, "] Start plotting spectrum", sep = "")
  chart.type <- tolower(trimws(chart.type))
  plot.name <- GenerateSampleName(obj)
  plot.title.suffix <- paste0("of Base Substitutions for \nSample \"",
                              plot.name, "\"")
  plot.out.name <- tolower(paste0(plot.name, ".spectrum.",
                                  chart.type, ".", class.name))
  
  if (chart.type == "num") {
    plot.data <- data.frame(obj$mut.base[, "count"], row.names = kMutBase)
    plot.title <- paste("Count", plot.title.suffix)
    cat(" (in Number) ...\n")
    
    ggplot(plot.data, aes(x = row.names(plot.data), y = plot.data)) +
      geom_bar(stat = "identity", fill = rev(kMutBaseColour), width = 0.4) +
      scale_x_discrete(limits = rev(row.names(plot.data))) +
      scale_y_continuous(labels = scales::comma) +
      coord_flip() +
      labs(title = plot.title, y = "Number of Mutations") +
      plotter$theme
  } else if (chart.type == "pct") {
    plot.data <- data.frame(obj$mut.base[, "percentage"], row.names = kMutBase)
    plot.title <- paste("Percentage", plot.title.suffix)
    cat(" (in Percentage) ...\n")
    
    ggplot(plot.data, aes(x = row.names(plot.data), y = plot.data)) +
      geom_bar(stat = "identity", fill = kMutBaseColour, width = 0.4) +
      scale_y_continuous(breaks = seq(0, max(plot.data), by = 0.05),
                         labels = scales::percent) +
      ggtitle(plot.title) +
      plotter$theme
  } else {
    stop("Chart type is not supported", call. = F)
  }
  
  ggsave(paste0(plotter$out.path, plot.out.name, ".", plotter$type))
  cat("[", class.name, "] Plotting completed in ", toupper(plotter$type), "\n",
      sep = "")
}
