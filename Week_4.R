# Title     : Exercises/tutorial BIF-30806
# Objective : TODO
# Created by: hidde
# Created on: 17/11/17

### DESeq2 R practical

library(DESeq2)
source("http://www.bioinformatics.nl/courses/BIF-30806/DEseq2Exercise.R")
internode_data=read.table(
  "http://www.bioinformatics.nl/courses/BIF-30806/maize_e3.table",
  row.names=1, header=TRUE, sep ="\t")

new_internode_data = internode_data[mx > 10,]
dim(new_internode_data)

count_data = new_internode_data[,1:6]
condition = factor(c("first","first","first","fourth","fourth","fourth"))
col_data = DataFrame(condition)

dds = DESeqDataSetFromMatrix(count_data, col_data, ~condition)
dds = estimateSizeFactors(dds)
sizeFactors(dds)

norm_versus_non_norm(dds, 1, 2, left = 2, right = 8)
rld = rlog(dds)
plot(density(assay(dds)[,1]), main="counts")
plot(density(assay(rld)[,1]), main="log counts")
dists = dist(t(assay(rld)))
plot(hclust(dists))

dds = estimateDispersions(dds)
plotDispEsts(dds)
dds = nbinomWaldTest(dds)
res = results(dds)

res$padj = ifelse(is.na(res$padj), 1, res$padj)
res$annotation = new_internode_data[,7]
write.table(res, col.names=NA, row.names=T, file ="internodes.tsv", sep ="\t")
plotMA(res, main="MA plot",ylim=c(-8,8),alpha=0.01)
