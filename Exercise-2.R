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


# Exercise to practice

https://rawgit.com/bioinformatics-core-shared-training/microarray-analysis/master/de-tutorial.html

https://rawgit.com/bioinformatics-core-shared-training/microarray-analysis/master/affymetrix.nb.html

library(affy)
setwd("C:/Users/Bio-user/Documents/GitHub/Class-Exercise-1/estrogen")
targetsFile <- "estrogen.txt"
pd <- read.AnnotatedDataFrame(targetsFile,header=TRUE,sep="",row.names=1)
pData(pd)


ER <- pData(pd)$estrogen
Time <- factor(pData(pd)$time.h)
design <- model.matrix(~ER+Time)
design

design2 <- model.matrix(~ER*Time)
design2

raw <-ReadAffy(celfile.path = "C:/Users/Bio-user/Documents/GitHub/Class-Exercise-1/estrogen", filenames=rownames(pData(pd)),phenoData = pd)
raw

boxplot(raw,col="red",las=2)

par(mfrow=c(2,1))
hist(log2(pm(raw[,1])),breaks=100,col="steelblue",main="PM",xlim=c(4,14))
hist(log2(mm(raw[,1])),breaks=100,col="steelblue",main="MM",xlim=c(4,14))


source("https://bioconductor.org/biocLite.R")
biocLite("affyPLM")

library(affyPLM)

plmset <- fitPLM(raw)
NUSE(plmset,las=2)

bad <- ReadAffy(celfile.path = "C:/Users/Bio-user/Documents/GitHub/Class-Exercise-1/estrogen",filenames="bad.cel")
image(bad)

par(mfrow=c(2,4))
image(raw[,1])

eset <- rma(raw)

library(limma)
fit1 <- lmFit(eset, design)
fit1 <- eBayes(fit1)
topTable(fit1, coef=2)


fit2 <- lmFit(eset, design2)
fit2 <- eBayes(fit2)
topTable(fit2, coef=2)

#Annotation of samples

library(GEOquery)
library(limma)
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33126/matrix/GSE33126_series_matrix.txt.gz"
filenm <- "GSE33126_series_matrix.txt.gz"
if(!file.exists(filenm)) download.file(url, destfile=filenm)
gse <- getGEO(filename=filenm)

# Download in Linux ##########################
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33126/matrix/GSE33126_series_matrix.txt.gz
gunzip GSE33126_series_matrix.txt.gz
Change permission 
chmod 777 All
To set user (owner) executable permission bit on:
chmod u+x file
To set group read / write permission bits:
chmod g+rw file
###################################

gse
head(exprs(gse))

exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)

pd <- pData(gse)
SampleGroup <- pd$source_name_ch1
library(genefilter)
destats <- rowttests(exprs(gse),fac=SampleGroup)
head(destats)

library(genefilter)
gse.expFilt <- varFilter(gse)
gse.expFilt

# Annotation

anno <- fData(gse.expFilt)
head(anno)[,1:5]

anno <- anno[,c("Symbol","Entrez_Gene_ID","Chromosome","Cytoband")]
fit2$genes <- anno
topTable(fit2)

volcanoplot(fit2)

testResults <- topTable(fit2, number=nrow(fit2))
testResults[which(testResults$Symbol == "OCIAD2"),]

mylist <- c("LOC441066","ARF3","FMNL3","CSF1R","XLKD1","TTRAP","DMWD","SYNPO2L","PILRB","LAMP3")
testResults[which(testResults$Symbol %in% mylist),]






































































