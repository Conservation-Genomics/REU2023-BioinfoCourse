####Script for Class DE
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

install.packages("statmod")
install.packages("corrplot")

#loading libraries
library("edgeR")
library("corrplot") 

#change working dir
setwd("~/Desktop/RClass/")

##Read in counts file
tomct <- read.delim("gene_count_matrix.csv", row.names=1, stringsAsFactors=FALSE, sep = ",")

#make groups
tom_group <- (c("im","im", "br", "br"))
tom_dgl <- DGEList(counts=tomct, group=tom_group)

#access components of DGE object
head(tom_dgl$counts)
dim(tom_dgl$counts)
head(tom_dgl$samples)

#store a copy of raw data in new variable
tom_dgl_raw <- tom_dgl

#Filter
keep <- rowSums(cpm(tom_dgl) > 10) >= 3  #keep genes with at least 10 cpm in at least 3  samples
tom_dgl <- tom_dgl[keep,]
table(keep)
dim(tom_dgl)

##Re-compute the library sizes:
tom_dgl$samples$lib.size <- colSums(tom_dgl$counts)

#check library sizes
png(filename = 'counts.png', width = 1500, height = 1500, units = 'px')
barplot(tom_dgl$samples$lib.size, main="Raw Counts", xlab = "sample", ylab = "counts")
dev.off()

# Get log2 counts per million
logcpm <- cpm(tom_dgl$counts,log=TRUE, prior.count = 1)

# Check distributions of samples using boxplots
# Let's add a blue horizontal line that corresponds to the median logCPM
#color by group
par(mfrow=c(1,1),oma=c(2,0,0,0)) ### fix colors Suzy!!!!
group.col <- c("green", "green", "orange", "orange")#[tom_group]
png(filename = 'logcpm.png', width = 1500, height = 1500, units = 'px')
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,col=group.col,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs\n(coloured by groups)",cex.main=0.8)
dev.off()

#Estimate normalization factors, By default, calcNormFactors uses the TMM method and the sample whose 75%-ile (of library-scale-scaled counts) is closest to the mean of 75%-iles as the reference.
tom_dgl = calcNormFactors(tom_dgl)
tom_dgl$samples

logcounts <- cpm(tom_dgl,log=TRUE)
par(mfrow=c(1,2))
plotMD(logcounts,column=1)
abline(h=0,col="grey")
plotMD(tom_dgl,column = 2)
abline(h=0,col="grey")
par(mfrow=c(1,1))

#estimate dispersion
tom_dgl <- estimateCommonDisp(tom_dgl, verbose=TRUE)
tom_dgl <- estimateTagwiseDisp(tom_dgl)
summary(tom_dgl$prior.df)
sqrt(tom_dgl$common.disp)  #The square root of the common dispersion gives the coefficient of variation of biological variation (BCV).
png(filename = 'bcv.png', width = 800, height = 800, units = 'px')
plotBCV(tom_dgl)
dev.off()

########QC############
png(filename = 'dendrogram.png', width = 800, height = 800, units = 'px')
tcounts <- t(as.table(as.matrix(tom_dgl)))
counts.dist = hclust(dist(tcounts))
plot(counts.dist)
dev.off()

#scatterplot of reps
#png(filename = 'pnp1_pnp2_scatter.png', width = 800, height = 800, units = 'px')
tom_dgl.cpm <-cpm(tom_dgl, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)
plot(tom_dgl.cpm[,"SRR404331_ch4.sort"], tom_dgl.cpm[,"SRR404333_ch4.sort"], log="xy")
#dev.off()

#Correlation Matrix
png(filename = 'corr_matrix.png', width = 1500, height = 1500, units = 'px')
corrplot(cor(tom_dgl.cpm), method="square", cl.lim=c(0,1), tl.col="black", addgrid.col="black", is.corr=FALSE, number.cex = 0.5)
dev.off()

#############DE#######################
#Set up comparisons
# Get all genes that are de in cfm to wt comparisons
tom_de <- exactTest(tom_dgl)

#get top 10 DE
topTags(tom_de,20)

#Select differentially expressed genes at a FDR of 5%
tom_de_0.5 <- decideTestsDGE(tom_de)
tom_de_0.5
summary(tom_de_0.5)

#Generate a DF with DE genes at FDR of 5%
isDE <- as.logical(tom_de_0.5)
tom_de_05.name <- rownames(tom_dgl)[isDE]
tom_de_05.table <- tom_de[tom_de_05.name,]
write.csv(tom_de_05.table, file="test.csv")

#Smear Plot
plotSmear(tom_de, de.tags = tom_de_05.name)
abline(h = c(-1,1), col = "blue")

##Install heatmap.2 packages
install.packages("gplots")
library(gplots)

##Plot heatmap
colors <- colorRampPalette(c("yellow", "orange", "red"))(n = 299)
png(filename = "logcpm.png", width = 700, height = 1000, units = "px", pointsize = 15)
  heatmap <- heatmap.2(as.matrix(logcpm), col=colors, na.rm=TRUE, labRow = NULL, dendrogram="row",
             Colv="NA", key=TRUE, symkey=FALSE, density.info="none", trace="none",
             margins=c(10,10), cexCol=1, symbreaks=FALSE)
dev.off()

