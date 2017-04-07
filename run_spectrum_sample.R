args <- commandArgs(trailingOnly = T)
#args <- c("Test.vcf.gz", "refseq.fasta", "pdf", "summary_out/")

# Constants
kMinArgLength <- 2

if (length(args) < kMinArgLength) {
  stop("Missing arguments", call. = F)
}

if (!exists("VariantSample", mode = "function")) {
  source("variant_sample.R")
}

cat("Loading VCF file: \"", basename(args[1]), "\" ...\n", sep = "")
raw.vcf <- tryCatch({
  readVcf(args[1], "hg19", ScanVcfParam(info = NA, samples = NA))
}, error = function(e) {
  stop("Failed to load this VCF file. Stop", call. = F)
})
cat("Loading completed\n")

if (!exists("XPlotter", mode = "function")) {
  source("plotter.R")
}

sample.name <- strsplit(basename(args[1]), split = "\\.")[[1]][1]
singlessm <- SingleSSM(raw.vcf, FaFile(args[2]), name = sample.name)

# Plot single sample spectrum in number
num.theme.extra <- kBarTheme + kBarTheme.num
if (length(args) > kMinArgLength &&
    tolower(args[kMinArgLength + 1]) %in% c("jpg", "jpeg")) {
  num.plotter <- JpgPlotter(plot.theme.extra = num.theme.extra)
} else {
  num.plotter <- PdfPlotter(plot.theme.extra = num.theme.extra)
}
PlotSpectrum(singlessm, plotter = num.plotter, chart.type = "num")

# Plot single sample spectrum in percentage
pct.theme.extra <- kBarTheme + kBarTheme.pct
if (length(args) > kMinArgLength &&
    tolower(args[kMinArgLength + 1]) %in% c("jpg", "jpeg")) {
  pct.plotter <- JpgPlotter(plot.theme.extra = pct.theme.extra)
} else {
  pct.plotter <- PdfPlotter(plot.theme.extra = pct.theme.extra)
}
PlotSpectrum(singlessm, plotter = pct.plotter)

if (length(args) > kMinArgLength + 1) {
  Summary(singlessm, out.path = args[kMinArgLength + 2])
}
