setwd("C:/workspace/ml/all_c_elegans")
library(e1071)
library(Matrix)
library(SparseM)

install.packages("SparseM",depen=TRUE) 

df <- read.table(file = 'NeuronalGeneCount', skipNul = T, header = TRUE, sep=",", quote = "")

df <- readMM(file='NeuronalGeneCount')

df <- read.matrix.csr(file = 'NeuronalGeneCount')

df <- load(file = 'NeuronalGeneCount')
?load

x <- load('NeuronalGeneCount')



x<-get(load('data/NeuronalGeneCount'))

write.csv(as.matrix(x),file="data/NeuronalGeneCount.csv",row.names=TRUE) # drops the rownames
row.names(x)
nrow(x)
ncol(x)
