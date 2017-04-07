Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
library(SomaticSignatures)

args <- args <- commandArgs(trailingOnly = T)
#args <- c("Test.vcf.gz", "refseq.fasta")

vcf <- readVcf(args[1], "hg19", ScanVcfParam(info = NA, samples = NA))
idx <- which(isSNV(vcf))

seq.ctx <- mutationContext(as(vcf[idx], "VRanges"), FaFile(args[2]), strand = T)
sampleNames(seq.ctx) <- "2"
motif.num <- t(motifMatrix(seq.ctx, group = "sampleNames", normalize = F))

write.table(motif.num, file = "check.tsv", quote = F, sep = "\t", col.names = NA)
