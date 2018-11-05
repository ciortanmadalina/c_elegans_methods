rm(list = ls())
.rs.restartR()

library(Matrix)
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)

df<-get(load('all_c_elegans/data/NeuronalGeneCount'))
rm(NeuronalGenCount)
df <- as.matrix(df)
anno <- colnames(df)
head(df[ , 1:3])

umi <- SingleCellExperiment(
  assays = list(counts = as.matrix(df)), 
  colData = anno
)

keep_feature <- rowSums(counts(umi) > 0) > 0
umi <- umi[keep_feature, ]
umi

