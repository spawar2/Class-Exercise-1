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

# Making MA plots 
# M: difference in average log intensities
# A: average log intensities

d <- rowMeans(mas5.e[,1:3]) - rowMeans(mas5.e[,4:6])
a <- rowMeans(mas5.e)
plot(a, d, ylim = c(-5, 5), main = "MAS 5.0 MA plot", xlab = "A", ylab = "M", pch = ".")
abline(h = c(-1, 1))

#Finding specific probesets in the plot

spikedin <- colnames(pData(spikein95))
spikedIndex <- match(spikedin,featureNames(mas5.eset))
points(a[spikedIndex], d[spikedIndex],pch = 19, col = "red")

# Trying RMA normalization

rma.eset <- rma(spikein95)
rma.e <- exprs(rma.eset)

# Making MA plots 

d <- rowMeans(rma.e[,1:3]) - rowMeans(rma.e[,4:6])
a <- rowMeans(rma.e)
plot(a, d, ylim = c(-2, 2), main = "RMA MA plot", xlab = "A", ylab = "M", pch = ".")
abline(h = c(-1, 1))

#Finding specific probesets in the plot

points(a[spikedIndex], d[spikedIndex],pch = 19, col = "red")

# Volcano Plots
source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")
library("genefilter")
pData(rma.eset) <- pData(mas5.eset)
tt <- rowMeans(rma.e)
lod <- tt
plot(d, lod, cex = 0.25, main = "Volcano plot for MA", xlim = c(-2, 2), xlab = "M", ylab = "A", yaxt = "n")
axis(2, at = seq(0, 3, by = 1), labels = 10^(-seq(0, 3, by = 1)))
points(d[spikedIndex], lod[spikedIndex], pch = 19, col = "red")
abline(h = 2, v = c(-1, 1))






































