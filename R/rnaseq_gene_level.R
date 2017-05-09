#---
#title: "rnaseq_gene_level.Rmd"
#author: "Yang Hu"
#date: "2017/5/5"
#output: html_document
#---
    ## Introduction
    
#In this section, we will focus on comparing the expression levels of genes across different samples. 
## Counting reads in genes

#In this section we will examine 75 samples from Waldenstrom Macroglobulinemia patients, 
#in a cohort called ???zhunter,??? from Harvard University. 

## Construct sample.table


setwd("/Users/yah2014/Dropbox/Public/Olivier/R/zhunter/Mutation");getwd() #set_up_enviroment

sample_zhunter.table <- read.csv("sample_table_zhunter.csv")
fileName <-paste0(sample_zhunter.table$BioSample, ".bam.count")
sample_zhunter.table <- data.frame(sampleName = sample_zhunter.table$sampleName,
                                   fileName = fileName,
                                   Condition= sample_zhunter.table$Condition,
                                   QC = sample_zhunter.table$QC,
                                   Label= sample_zhunter.table$Label,
                                   BioSample = sample_zhunter.table$BioSample)
head(sample_zhunter.table)


## Creating a DESeqDataSet object

library(DESeq2)

directory <- "/Users/yah2014/Dropbox/Public/Olivier/R/ALL_COUNTS/zhunter_Counts" #dir for HTSeq Counts
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sample_zhunter.table,
                                       directory = directory,
                                       design= ~ Condition)
ddsHTSeq


#if Error in `colnames<-`(`*tmp*`, value = 1:75) : 
#    attempt to set 'colnames' on an object with less than two dimensions
#clean the first row with the empty first column and "0" in the second column in HTSeq.count files
#use terminal, run below command lines to remove first line for all files:
#    cd /Users/yah2014/Dropbox/Public/Olivier/R/ALL_COUNTS/zhunter_Counts
#for file in $(ls);do echo "$(tail -n +2 $file)" > $file;done

### Normalization for sequencing depth
####Pre-filtering

#removing rows in which there are no reads or nearly no reads

ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ]


#Note on factor levels


ddsHTSeq$Condition <- factor(ddsHTSeq$Condition, levels=c("Cancer","Healthy","Unclear"))
ddsHTSeq$Label <- factor(ddsHTSeq$Label, levels=c("NWM","rnaWT","WM","MWCL","BCWM"))


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
        
plotPCA(vsd, intgroup="Condition")

    
#We can make this plot even nicer using custom code from the *ggplot2* library:
        
library(ggplot2)
head((data <- plotPCA(vsd, intgroup=c("Condition","Label"), returnData=TRUE)))
    (percentVar <- 100*round(attr(data, "percentVar"),2))

    
makeLab <- function(x,pc) paste0("PC",pc,": ",x,"% variance")
ggplot(data, aes(PC1,PC2,col=Condition,shape=Label)) + geom_point() +
        xlab(makeLab(percentVar[1],1)) + ylab(makeLab(percentVar[2],2))

    
#In addition, we can plot a hierarchical clustering based on Euclidean distance matrix:
        

plot(hclust(dist(t(log.norm.counts))), labels=colData(dds)$Condition,cex = 0.75)
plot(hclust(dist(t(assay(vsd)))), labels=colData(vsd)$Condition,cex = 0.75)

## Differential gene expression
### Modeling raw counts with normalization
    
### Experimental design and running DESeq2
dds <- DESeqDataSet(ddsHTSeq, design= ~ Condition)

#First, we setup the `design` of the experiment, so that differences will be considered across time and protocol variables.
#We can read and if necessary reset the design using the following code.
        

design(dds)
design(dds) <- ~Condition

#The last variable in the design is used by default for building results tables 
#(although arguments to `results` can be used to customize the results table), 
#and we make sure the "control" or "untreated" level is the first level, 
#such that log fold changes will be treated over control, and not control over treated.
        
        
levels(dds$Condition)
dds$Condition <- relevel(dds$Condition, "Unclear")
dds$Condition <- relevel(dds$Condition, "Healthy")
levels(dds$Condition)

        
#The following line runs the *DESeq2* model. 
#After this step, we can build a results table, 
#which by default will compare the levels in the last variable in the design, 
#so the *Condition* treatment in our case:
            
dds <- DESeq(dds)
res <- results(dds)

### Examining results tables
        
head(res)
table(res$padj < 0.1)

#A summary of the results can be generated:
            
summary(res)

        
#For testing at a different threshold, we provide the `alpha` to *results*,
#so that the mean filtering is optimal for our new FDR threshold.
        
res2 <- results(dds, alpha=0.05)
table(res2$padj < 0.05)

#Exporting results to CSV files
        
resO5rdered <- res2[order(res2$padj),]
write.csv(as.data.frame(resO5rdered), 
file="condition_treated_results.csv")

### Visualizing results
        
#The MA-plot provides a global view of the differential genes, 
#with the log2 fold change on the y-axis over the mean of normalized counts:
            
plotMA(res, ylim=c(-8,4))

#We can also test against a different null hypothesis. 
#For example, to test for genes which have fold change more than doubling or less than halving:
            
res.thr <- results(dds, lfcThreshold=1)
plotMA(res.thr, ylim=c(-8,4))

        
        
#A p-value histogram:
            
hist(res$pvalue[res$baseMean > 1], 
col="grey", border="white", xlab="", ylab="", main="")

#A sorted results table:
            
resSort <- res[order(res$padj),]
head(resSort)

#Examine the counts for the top gene, sorting by p-value:
 

#plotCounts(dds, gene=which.min(res$padj), intgroup="Condition")

#Make normalized counts plots for the top 9 genes:
            
par(mfrow=c(3,3))
for (i in 1:9)  plotCounts(dds, order(res$padj)[i], intgroup="Condition")

#A more sophisticated plot of counts:
            
library(ggplot2)
data <- plotCounts(dds, gene=which.min(res$padj), intgroup=c("Condition","Label"), returnData=TRUE)
ggplot(data, aes(x=Condition, y=count, col=Label)) +
        geom_point(position=position_jitter(width=.1,height=0)) +
        scale_y_log10()

#Connecting by lines shows the differences which are actually being tested by *results* given that our design includes `cell + Condition`
        
par(mfrow=c(1,1))
ggplot(data, aes(x=Condition, y=count, col=Label, group=Label)) +
       geom_point() + geom_line() + scale_y_log10() 

#A heatmap of the top genes:
            
library(pheatmap)
topgenes <- head(rownames(resSort),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Condition","Label")])
pheatmap(mat, annotation_col=df, fontsize_col=6)

        
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
plot(svseq$sv[,1], svseq$sv[,2], col=dds$Label,pch=16)
legend("topright",pch = 16, col=1:5,levels(dds$Label))
#text(svseq$sv[,1], svseq$sv[,2], 1:ncol(dds), pos=1)

#Do the surrogate variables capture the health condition?
        
plot(svseq$sv[,1], svseq$sv[,2], col=dds$Condition,pch=16)
legend("topright",pch = 16, col=1:3,levels(dds$Condition))
#text(svseq$sv[,1], svseq$sv[,2], 1:ncol(dds), pos=1)

        
#Using the surrogate variables in a *DESeq2* analysis:

#dds.sva <- dds
#dds.sva$SV1 <- svseq$sv[,1]
#dds.sva$SV2 <- svseq$sv[,2]
#design(dds.sva) <- ~ SV1 + SV2 + Condition
#dds.sva <- DESeq(dds.sva)

        
## Session info
    
sessionInfo()
