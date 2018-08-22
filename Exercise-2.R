# R Exercise 2: Analysis of Microarray data, Shrikant Pawar 08/22/2018

source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("SpikeInSubset")
library(SpikeInSubset)
data(spikein95)
image(spikein95[, 3])



