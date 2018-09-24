# R Exercise 2: Analysis of Microarray data, Shrikant Pawar 08/22/2018
# function to separate two points, Intercept as 0??

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
# Additive Effect, no difference between 10 and 48
##   (Intercept) ER 48h
## 1           1     0    0
## 2           1     0    1
# For any gene1 if it finds P value <0.05, it will be significant

design2 <- model.matrix(~ER*Time)
design2
# Interactions, accounts for differences in time
##   (Intercept)    ER   48h     ER:48h
## 1           1     0    0          0
## 2           1     0    1          0
# If gene1 P value <0.05 with 10hr, checks for same with 48hr

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

#Annotation of the FIT samples
setwd("C:/Users/Bio-user/Documents/GitHub/Class-Exercise-1/estrogen")
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

##########Not used below #########################

testResults <- topTable(fit2, number=nrow(fit2))
testResults[which(testResults$Symbol == "OCIAD2"),]

mylist <- c("LOC441066","ARF3","FMNL3","CSF1R","XLKD1","TTRAP","DMWD","SYNPO2L","PILRB","LAMP3")
testResults[which(testResults$Symbol %in% mylist),]

#####################################################

In Class Assignment:

# Load the libraries GEOquery, affy and limma.

library(GEOquery)
library(limma)
library(affy)

# Use url FTP link to download sample dataset GSE1000_series_matrix.txt.gz
with function getGEO and assign it to an object gse.

url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/matrix/GSE1000_series_matrix.txt.gz"
filenm <- "GSE1000_series_matrix.txt.gz"
if(!file.exists(filenm)) download.file(url, destfile=filenm)
gse <- getGEO(filename=filenm)

# Check the initial few rows of object gse with a function head.

head(exprs(gse))

# Take a log to base 2 for all expression levels

gse_log <- log2(exprs(gse))

# First 5 columns in gse are treatment, and next 5 are control.Take a average
of first five and put in object named treatment, take average of last 5 columns
and put in object named control. Then take a simple fold change to calculate
fold difference in your expression levels and put it in object named fold.
Use gse_log object to divide treatment and control columns use function 
rowMeans() to calcuate the means of rows.

treatment <- rowMeans(gse_log[,1:5])
control <- rowMeans(gse_log[,6:10])
fold <- treatment/control

# Create box plot for treatment object in green and control object in red color
on one panel for the calculated row means!!
par(mfrow=c(2,1))
boxplot(treatment, col="green")
boxplot(control, col="red")

# Export the fold object to a csv file.

write.csv(fold, file = "fold.csv",row.names=FALSE)























 





















































































