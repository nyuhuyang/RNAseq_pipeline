#---
#title: "rnaseq_gene_level.Rmd"
#author: "Yang Hu"
#date: "2017/5/5"
#output: html_document
#---
## Introduction

#In this section, we will focus on comparing the expression levels of genes across different samples. 
## Counting reads in genes

#In this section we will examine 80 samples from Waldenstrom Macroglobulinemia patients, 
#in a cohort called zhunter from France. 


#detect OS and set enviroment
if (Sys.info()[['sysname']]=="Darwin"){
    setwd("/Users/yah2014/Dropbox/Public/Olivier/R/zhunter/Mutation");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
    setwd("C:/Users/User/Dropbox/Public/Olivier/R/zhunter/Mutation");getwd();list.files()}

# Construct sample.table
sample_zhunter.table <- read.csv("sample_table_zhunter.csv")
fileName <-paste0(sample_zhunter.table$Sample_ID, ".bam.count")
sample_zhunter.table <- data.frame(sampleName = sample_zhunter.table$Sample_name,
                                 fileName = fileName,
                                 Condition= sample_zhunter.table$Condition,
                                 SPI1=sample_zhunter.table$SPI1,
                                 Sample_ID = sample_zhunter.table$Sample_ID)
head(sample_zhunter.table)


## Creating a DESeqDataSet object

library(DESeq2)
#detect OS and asign dir for HTSeq Counts
if (Sys.info()[['sysname']]=="Darwin"){
    directory <- "/Users/yah2014/Dropbox/Public/Olivier/R/ALL_COUNTS/zhunter_Counts"}
if (Sys.info()[['sysname']]=="Windows"){
    directory <- "C:/Users/User/Dropbox/Public/Olivier/R/ALL_COUNTS/zhunter_Counts"}

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sample_zhunter.table,
                                       directory = directory,
                                       design= ~ Condition) # Don't +SPI1
ddsHTSeq
#if Error in `colnames<-`(`*tmp*`, value = ....) : 
#    attempt to set 'colnames' on an object with less than two dimensions
#Way to solve:
#clean the first row with the empty first column and "0" in the second column in HTSeq.count files
#use terminal, run below command lines to remove first line for all files:
#    cd /Users/yah2014/Dropbox/Public/Olivier/R/ALL_COUNTS/zhunter_Counts
#for file in $(ls);do echo "$(tail -n +2 $file)" > $file;done

### Normalization for sequencing depth
####Pre-filtering
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) !=0, ] #removing rows with 0 reads


#Note on factor levels


ddsHTSeq$Condition <- factor(ddsHTSeq$Condition, levels=c("Healthy","Cell line","Unclear","Cancer"))
ddsHTSeq$SPI1 <- factor(ddsHTSeq$SPI1, levels=c("um","m")) 

#The following estimates size factors to account for differences in sequencing depth, 
#and is only necessary to make the `log.norm.counts` object below.


dds <- estimateSizeFactors(ddsHTSeq)
head(sizeFactors(dds))
head(colSums(counts(dds)))
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))


#Size factors are calculated by the median ratio of samples to
#a pseudo-sample (the geometric mean of all samples). 
#In other words, for each sample, we take the exponent of the median of the log ratios in this histogram.


loggeomeans <- rowMeans(log(counts(dds)))
hist(log(counts(dds)[,1]) - loggeomeans, 
     col="grey", main="", xlab="", breaks=40)


#The size factor for the first sample:

exp(median((log(counts(dds)[,1]) - loggeomeans)[is.finite(loggeomeans)]))
sizeFactors(dds)[1]

#Make a matrix of log normalized counts (plus a pseudocount):

log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)

#Examine the log normalized counts (plus a pseudocount).


rs <- rowSums(counts(dds))
boxplot(log.norm.counts[rs > 0,]) # normalized


#Make a scatterplot of log normalized counts against each other.
#Note the fanning out of the points in the lower left corner, for points less than $2^5 = 32$.


plot(log.norm.counts[,1:2], cex=.1)


### Stabilizing count variance

#Transformation for stabilizing variance in the *DESeq2* package is `varianceStabilizingTransformation`.
#These two tranformations are similar, the *rlog* might perform a bit better when the size factors vary widely,
#and the *varianceStabilizingTransformation* is much faster when there are many samples.

vsd <- varianceStabilizingTransformation(dds)
plot(assay(vsd)[,1:2], cex=.1)


# We can examine the standard deviation of rows over the mean for the *vsd*.
#Note that the genes with high variance for the *log* come from the genes with lowest mean.
#If these genes were included in a distance calculation, 
#the high variance at the low count range might overwhelm the signal at the higher count range.

library(vsn)
meanSdPlot(assay(vsd), ranks=FALSE)


#The principal components (PCA) plot is a useful diagnostic for examining relationships between samples:

#Using the VST:

plotPCA(vsd, intgroup="SPI1")



#We can make this plot even nicer using custom code from the *ggplot2* library:

library(ggplot2)
head((data <- plotPCA(vsd, intgroup=c("Condition","SPI1"), returnData=TRUE)))
(percentVar <- 100*round(attr(data, "percentVar"),2))


makeLab <- function(x,pc) paste0("PC",pc,": ",x,"% variance")
g<-ggplot(data, aes(PC1,PC2,col=Condition,shape=SPI1))
g+ geom_point(aes(size = SPI1)) +
    xlab(makeLab(percentVar[1],1))+ 
    ylab(makeLab(percentVar[2],2))+
    labs(title = "PCA for Health Condition and SPI1 mutation type")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),text = element_text(size=32))

#In addition, we can plot a hierarchical clustering based on Euclidean distance matrix:

plot(hclust(dist(t(log.norm.counts))), SPI1s=colData(dds)$Condition,cex = 0.75)
plot(hclust(dist(t(assay(vsd)))), SPI1s=colData(vsd)$Condition,cex = 0.75)

## Differential gene expression
### Modeling raw counts with normalization

### Experimental design and running DESeq2
dds <- DESeqDataSet(ddsHTSeq, design= ~ Condition+SPI1)

#First, we setup the `design` of the experiment, so that differences will be considered across time and protocol variables.
#We can read and if necessary reset the design using the following code.


design(dds)
design(dds) <- ~Condition+SPI1

#The last variable in the design is used by default for building results tables 
#(although arguments to `results` can be used to customize the results table), 
#and we make sure the "control" or "untreated" level is the first level, 
#such that log fold changes will be treated over control, and not control over treated.


levels(dds$Condition)
levels(dds$SPI1)

#The following line runs the *DESeq2* model. 
#After this step, we can build a results table, 
#which by default will compare the levels in the last variable in the design, 
#so the *Condition* treatment in our case:

dds <- DESeq(dds)
res_cancer <- results(dds,contrast = c("Condition","Cancer","Healthy"))
res_SPI1 <- results(dds,contrast = c("SPI1","m","um"))
### Examining results tables

head(res_cancer)
table(res_cancer$padj < 0.1)
head(res_SPI1)
table(res_SPI1$padj < 0.1)
#A summary of the results can be generated:

summary(res_cancer)


#For testing at a different threshold, we provide the `alpha` to *results*,
#so that the mean filtering is optimal for our new FDR threshold.

res2 <- results(dds, alpha=0.5,lfcThreshold=1)
table(res2$padj < 0.05)

#Exporting results to CSV files

resO5rdered <- res2[order(res2$padj),]
resO5rdered <- resO5rdered[complete.cases(resO5rdered),] #remove NA
resE8 <- resO5rdered[resO5rdered$padj<0.05,]
write.csv(as.data.frame(resE8),file="condition_treated_results_Condition_zhunter.csv")
#RUN GESA on https://david.ncifcrf.gov/ by uploading gene names


### Visualizing results

#The MA-plot provides a global view of the differential genes, 
#with the log2 fold change on the y-axis over the mean of normalized counts:

plotMA(res_cancer, ylim=c(-5,5))

#We can also test against a different null hypothesis. 
#For example, to test for genes which have fold change more than doubling or less than halving:

res.thr <- results(dds, lfcThreshold=1)
plotMA(res.thr, ylim=c(-5,5))



#A p-value histogram:

hist(res_cancer$pvalue[res_cancer$baseMean > 1], 
     col="grey", border="white", xlab="", ylab="", main="")


#Examine the counts for the top gene, sorting by p-value:

par(mfrow=c(1,1))
par(oma=c(2,2,2,2))
plotCounts(dds, gene=which.min(res_cancer$padj), intgroup="Condition")
plotCounts(dds, which(rownames(dds)=="SPI1"), intgroup="SPI1", 
           main="SPI1 expression level in unmutated vs mutated group",
           cex.main = 3,cex=1.5,cex.axis=1.5)

d <- plotCounts(dds, which(rownames(dds)=="SPI1"), intgroup="SPI1",returnData=TRUE)
rownames(d[which(d$count<200),])

#Make normalized counts plots for the top 9 genes:
par(mfrow=c(3,3))
for (i in 1:9)  plotCounts(dds, order(res_cancer$padj)[i], intgroup="SPI1")

#A more sophisticated plot of counts:

library(ggplot2)
data <- plotCounts(dds, gene=which(rownames(res_cancer)=="SPI1"), intgroup=c("Condition","SPI1"), returnData=TRUE)
ggplot(data, aes(x=Condition, y=count, col=SPI1,size = SPI1))+
    geom_point(position=position_jitter(width=.1,height=0))+
    labs(title = "SPI1 expression level in Healthy vs Cancer group")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),text = element_text(size=40)) + 
    scale_y_log10()

ggplot(data, aes(x=SPI1, y=count, col=Condition,size = SPI1))+
    geom_point(position=position_jitter(width=.1,height=0))+
    labs(title = "SPI1 expression level in unmutated vs mutated group")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),text = element_text(size=38)) + 
    scale_y_log10()

#A sorted results table:

res_cancer_Sort <- res_cancer[order(res_cancer$padj),]
head(res_cancer_Sort)

res_SPI1_Sort <- res_SPI1[order(res_SPI1$padj),]
head(res_SPI1_Sort)
#A heatmap of the top genes:

library(pheatmap)
par(mfrow=c(1,1))
topgenes <- head(rownames(res_cancer_Sort),100)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Condition","SPI1")])
pheatmap(mat, annotation_col=df,
         fontsize=13,fontsize_row=7,fontsize_col=13,cex=1,
         main ="Top 100 differentially expressed genes between Healthy and Cancer")

topgenes <- head(rownames(res_SPI1_Sort),100)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Condition","SPI1")])
pheatmap(mat, annotation_col=df, 
         fontsize=13,fontsize_row=7,fontsize_col=13,cex=1,
         main ="Top 100 differentially expressed genes between mutated SPI1 and unmutated SPI1")
### Getting alternate annotations


#We can then check the annotation of these highly significant genes:


#library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
#anno <- select(org.Hs.eg.db, keys=topgenes,
#               columns=c("SYMBOL","GENENAME"), 
#               keytype="ENSEMBL")
#anno[match(topgenes, anno$ENSEMBL),]
#for Bioconductor >= 3.1, easier to use mapIds() function


### Looking up different results tables

#The `contrast` argument allows users to specify what results table should be built.
#See the help and examples in `?results` for more details:

#results(dds, contrast=c("cell","N61311","N052611"))


### Surrogate variable analysis for RNA-seq

#If we suppose that we didn't know about the different cell-lines in the experiment,
#but noticed some structure in the counts, 
#we could use surrograte variable analysis (SVA) to detect this hidden structure
#(see PH525x Course 3 for details on the algorithm).

library(sva)
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ Condition, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)


# Do the surrogate variables capture the cell difference?

par(mfrow=c(1,1))
par(oma=c(3,3,3,3))
par(mar=c(5,2,2,2))
plot(svseq$sv[,1], svseq$sv[,2], col=dds$SPI1,pch=16,cex=2.5,cex.axis=3,
     main="Batch effect analysis(sva) for zhunter cohort",cex.main = 3,
     sub="Mutated SPI1(red) clustered into separated batches",cex.sub = 3,cex.lab=0.01)
legend("topleft",pch = 16, col=1:5,cex=3,levels(dds$SPI1))
text(svseq$sv[,1], svseq$sv[,2], colnames(dds), pos=4)

#Do the surrogate variables capture the health condition?

plot(svseq$sv[,1], svseq$sv[,2], col=dds$Condition,pch=16,main="Surrograte Variable Analysis (SVA)")
legend("topleft",pch = 16, col=1:3,levels(dds$Condition))
text(svseq$sv[,1], svseq$sv[,2], colnames(dds), pos=2)


#Using the surrogate variables in a *DESeq2* analysis:

#dds.sva <- dds
#dds.sva$SV1 <- svseq$sv[,1]
#dds.sva$SV2 <- svseq$sv[,2]
#design(dds.sva) <- ~ SV1 + SV2 + Condition
#dds.sva <- DESeq(dds.sva)


## Session info

sessionInfo()

