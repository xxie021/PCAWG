library(VariantAnnotation)
library(SomaticSignatures)

# Constants
kCategoryTypeI <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
kColourTypeI <- c("deepskyblue", "black", "tomato",
                  "gray", "yellowgreen", "pink")

# Inheritance
# Top: VariantSample
# 2nd Level: SsmSample, SimSample
# 3rd Level: SsmTypeI > SsmSample
VariantSample <- function(vcf, ref, id = "", name = "") {
  if (!(class(vcf) %in% c("CollapsedVCF", "ExpandedVCF"))) {
    stop("Invalid VCF data found", call. = F)
  }
  
  if (class(ref) != "FaFile") {
    stop("Invalid reference sequence found", call. = F)
  }
  
  tryCatch({
    countFa(ref)
  }, error = function(e) {
    stop("Invalid reference sequence found", call. = F)
  })
  
  #if (class(vcf) == "CollapsedVCF") {
  #  vcf = expand(vcf)
  #}
  
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

SsmSample <- function(vcf, ref, id = "", name = "") {
  me <- VariantSample(vcf, ref, id, name)
  idx <- which(isSubstitution(me$vcf))
  
  seq.ctx <- mutationContext(as(vcf[idx], "VRanges"), ref)
  sampleNames(seq.ctx) <- id
  #motif.pct <- motifMatrix(seq.ctx, group = "sampleNames")
  motif.num <- motifMatrix(seq.ctx, group = "sampleNames", normalize = F)
  #total.base.mutation <- sum(motif.num)
  #set.ref <- ref(me$vcf)[idx]
  #set.alt <- alt(me$vcf)[idx]
  
  me[["idx"]] <- idx
  me[["total.mut"]] <- sum(motif.num)
  #me[["seq.ctx"]] <- seq.ctx
  me[["c2a.ctx"]] <- CountBaseMutationWithContext(motif.num, "c", "a")
  me[["c2a.ctx"]] <- CountBaseMutationWithContext(motif.num, "c", "a")
  me[["c2a.ctx"]] <- CountBaseMutationWithContext(motif.num, "c", "a")
  me[["c2a.ctx"]] <- CountBaseMutationWithContext(motif.num, "c", "a")
  me[["c2a.ctx"]] <- CountBaseMutationWithContext(motif.num, "c", "a")
  me[["c2a.ctx"]] <- CountBaseMutationWithContext(motif.num, "c", "a")
  #me[["set.ref"]] <- set.ref
  #me[["set.alt"]] <- set.alt
  
  class(me) <- append(class(me), "SsmSample")
  return(me)
}

SsmTypeI <- function(vcf, id = "", name = "") {
  me <- SsmSample(vcf, id, name)
  
  me[["idx.c2a"]] <- GetIndexOfSsmTypeI(me, "c", "a")
  me[["idx.c2g"]] <- GetIndexOfSsmTypeI(me, "c", "g")
  me[["idx.c2t"]] <- GetIndexOfSsmTypeI(me, "c", "t")
  me[["idx.t2a"]] <- GetIndexOfSsmTypeI(me, "t", "a")
  me[["idx.t2c"]] <- GetIndexOfSsmTypeI(me, "t", "c")
  me[["idx.t2g"]] <- GetIndexOfSsmTypeI(me, "t", "g")
  
  idx.count <- c(length(me$idx.c2a), length(me$idx.c2g), length(me$idx.c2t),
                 length(me$idx.t2a), length(me$idx.t2c), length(me$idx.t2g))
  me[["out.data"]] <- data.frame(count = idx.count,
                                 percentage = sapply(idx.count, function(x) {
                                   x / sum(idx.count)
                                 }))
  
  class(me) <- append(class(me), "SsmTypeI")
  return(me)
}

CountBaseMutationWithContext <- function(obj, nt.ref, nt.alt) {
  UseMethod("CountBaseMutationWithContext", obj)
}

CountBaseMutationWithContext.matrix <- function(obj, nt.ref, nt.alt) {
  mut.base <- toupper(paste0(nt.ref, nt.alt))
  mut.matched.ctx <- obj[which(sapply(strsplit(row.names(obj), " "),
                                      function(x) x[1]) == mut.base), ]
  
  return(mut.matched.ctx)
}

GetIndexOfSsmTypeI <- function(obj, nt.ref, nt.alt,
                               complementarity = TRUE,
                               sort = FALSE) {
  UseMethod("GetIndexOfSsmTypeI", obj)
}

GetIndexOfSsmTypeI.SsmSample <- function(obj, nt.ref, nt.alt,
                                         complementarity = TRUE,
                                         sort = FALSE) {
  if (is.character(nt.ref)) {
    nt.ref <- DNAString(nt.ref)
  }
  if (is.character(nt.alt)) {
    nt.alt <- DNAString(nt.alt)
  }
  
  idx.ssm <- intersect(which(obj$set.ref == nt.ref),
                       which(obj$set.alt == nt.alt))
  if (complementarity) {
    idx.ssm.comp <- intersect(which(obj$set.ref == complement(nt.ref)),
                              which(obj$set.alt == complement(nt.alt)))
    idx.ssm <- c(idx.ssm, idx.ssm.comp)
  }
  
  if (sort) {
    idx.ssm <- sort(idx.ssm)
  }
  
  return(idx.ssm)
}

GenFileNamePrefix <- function(obj, full = FALSE) {
  UseMethod("GenFileNamePrefix", obj)
}

GenFileNamePrefix.VariantSample <- function(obj, full = FALSE) {
  if (full && obj$id != "" && obj$name != "") {
    return(paste0(obj$id, ".", obj$name))
  }
  
  return(ifelse(obj$id == "", obj$name, obj$id))
}

Summary <- function(obj, smry.out.path = NA) {
  if (!is.na(smry.out.path) && !dir.exists(smry.out.path)) {
    dir.create(smry.out.path)
  }
  UseMethod("Summary", obj)
}

Summary.SsmTypeI <- function(obj, smry.out.path = NA) {
  smry.data <- obj$out.data
  smry.data$percentage <- round(smry.data$percentage * 100, digits = 2)
  row.names(smry.data) <- c("C:G > A:T", "C:G > G:C", "C:G > T:A",
                            "T:A > A:T", "T:A > C:G", "T:A > G:C")
  
  cat("Sample ID:\t", obj$id, "\n")
  cat("Sample name:\t", obj$name, "\n")
  print(smry.data)
  cat("\n")
  
  if (!is.na(smry.out.path)) {
    write.table(smry.data, file = paste0(smry.out.path,
                                         tolower(GenFileNamePrefix(obj)),
                                         ".summary.tsv"),
                quote = F, sep = "\t", col.names = NA)
  }
}

PlotSignature.SsmTypeI <- function(obj, plotter, chart.type = "pct") {
  cat("Start plotting single sample SSM Type I signature")
  
  plot.name <- GenFileNamePrefix(obj)
  plot.title <- paste0("Type I Signature of\nSample \"", plot.name, "\"")
  plot.out.name <- paste0(tolower(plot.name), ".sig.ssm.typeI")
  chart.type <- tolower(trimws(chart.type))
  
  if (chart.type == "num") {
    plot.data <- data.frame(obj$out.data$count, row.names = kCategoryTypeI)
    plot.out.name <- paste0(plot.out.name, ".num")
    cat(" (in Number) ...\n")
    
    ggplot(plot.data, aes(x = row.names(plot.data), y = plot.data)) +
      geom_bar(stat = "identity", fill = rev(kColourTypeI), width = 0.4) +
      scale_x_discrete(limits = rev(row.names(plot.data))) +
      scale_y_continuous(labels = scales::comma) +
      coord_flip() +
      labs(title = plot.title, y = "Number of Mutations") +
      plotter$plot.theme + plotter$axis.theme
  } else if (chart.type == "pct") {
    plot.data <- data.frame(obj$out.data$percentage, row.names = kCategoryTypeI)
    plot.out.name <- paste0(plot.out.name, ".pct")
    cat(" (in Percentage) ...\n")
    
    ggplot(plot.data, aes(x = row.names(plot.data), y = plot.data)) +
      geom_bar(stat = "identity", fill = kColourTypeI, width = 0.4) +
      scale_y_continuous(breaks = seq(0, max(plot.data), by = 0.05),
                         labels = scales::percent) +
      ggtitle(plot.title) +
      plotter$plot.theme + plotter$axis.theme
  } else {
    stop("Chart type is not supported", call. = F)
  }
  
  ggsave(paste0(plotter$out.path, plot.out.name, ".", plotter$type))
  
  cat("Plotting completed in", toupper(plotter$type) ,"\n")
}
