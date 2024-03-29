args <- commandArgs(trailingOnly = TRUE)

#.libPaths(args[3])
#library(dplyr)
library(stringr)
#library(tidyr)


mosdepth = read.table(args[1], header=T)
mt_mean_coverage = mosdepth$mean[mosdepth$chrom=='MT']
wgs_mean_coverage= mosdepth$mean[mosdepth$chrom=='total']
haplocheck = read.table(args[2], header=T)
sample = unlist(strsplit(as.character(haplocheck$Sample), "_"))[1]
contamination = haplocheck$Contamination.Level
if (contamination == 'ND') {
	contamination = 0
} 


combined_table = cbind(sample, sample, contamination, wgs_mean_coverage, mt_mean_coverage)
write.table(combined_table, "conta_cov.tsv", row.names=F, col.names=F, quote=F, sep="\t")
