# R Exercise 2: Analysis of Microarray data, Shrikant Pawar 08/22/2018

source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("SpikeInSubset")
library(SpikeInSubset)
data(spikein95)
image(spikein95)
ids <- geneNames(spikein95)
ids[1:10]
mas5.eset <- mas5(spikein95)
mas5.e <- log2(exprs(mas5.eset))
boxplot(spikein95)
x11()
boxplot(mas5.e, col = 2:5)
density1 <- density(mas5.e[, 1])
plot(density1, main = "MAS5 expression measure distributions")
density2 <- density(mas5.e[, 2])
lines(density2, col = "red")
density3 <- density(mas5.e[, 3])
lines(density3, col = "blue")









