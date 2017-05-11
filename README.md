----------------------------------------------<br />
RNA-Seq Analysis Pipeline: detecting mutations in RNA-Seq samples, exam genotype vs phenotype relationship<br />
Author: Yang Hu<br />
Department of Physiology and Biophysics<br />
Weill Cornell Medicine<br />
Email: Hu.Yang@nyu.edu<br />
----------------------------------------------<br />

## I. Introduction

This RNA-Seq Analysis pipeline has three sections:

1) Alignment and Counts<br />
Including STAR alignment, Samtools sort, HTSeq Count, and Cufflinks. Linux bash shell script is here:https://github.com/nyuhuyang/RNA-Seq-Analysis/blob/master/Shell%20script/Alignment_and_Counts.sh

2) Call Somatic mutations and annotate vcf<br />
Including Mpileup, Varscan, and SnpEff. Linux bash shell script is here:https://github.com/nyuhuyang/RNA-Seq-Analysis/blob/master/Shell%20script/CallVar_and_Annotation.sh

3) Summarize FPKM, QC and cluster<br />
Including mutation data clean, FPKM data clean, QC, hclust, and mutation/expression comparison.
R code is here:https://github.com/nyuhuyang/RNA-Seq-Analysis/blob/master/R/Analyzing%20somatic%20mutations%20in%20RNA-seq%20data.R
R Markdown demostration is here: https://github.com/nyuhuyang/RNA-Seq-Analysis/blob/master/vignettes/Analyzing%20somatic%20mutations%20in%20RNA-seq%20data.Rmd

4) DESeq for gene level comparision<br />
R code is here: https://github.com/nyuhuyang/RNA-Seq-Analysis/blob/master/R/Analyzing%20somatic%20mutations%20in%20RNA-seq%20data.R
R Markdown demostration is here: https://github.com/nyuhuyang/RNA-Seq-Analysis/blob/master/vignettes/rnaseq_gene_level.Rmd


## II. Work flow

!![plot of chunk Flow_work](vignettes/Flow_work.png)


### 1). Alignment_and_Counts.sh

  Input Materials:<br />
  -i RNA-seq fastq files<br />
  
  Output Materials:<br />
  -o Aligned, sort, and indexed Bam files<br />
  -o HTSeq Counts<br />
  -o FPKM<br />
  
  This bash shell is written for analyzing multiple projects.
  
  There are total 155 Waldenström Macroglobulinemia (more than one terabyte) RNA-seq data.
  Those Waldenström Macroglobulinemia RNA-seq data are from two separate projects.<br />
  Project zhunter: 75 RNA-seq data from Harvard.<br />
  Project WTCHG: 80 RNA-seq data from France.<br />
  
  The STAR-HTSeq.sh is project specific. The current script is for project zhunter only.
  The PROJECT_NAME can be changed to fit different projects.
  
###  2). CallVar_and_Annotation.sh

  Input Materials:<br />
  -i Aligned, sorted, and indexed Bam files<br />
  
  Output Materials:<br />
  -o mpileup file<br />
  -o mutations.vcf<br />
  -o dbsnp_annotated.vcf<br />
  -o cosmic_dbsnp.vcf<br />
  -o Annotated.eff.vcf<br />
  -o snpEff_genes.txt<br />
  -o snpEff_summary.html<br />
  <br />
  
### 2-1). Between section 2 and 3, run following script to multiple Annotated.eff.vcf files with bcftools:
  #merge multiple Annotated.eff.vcf files with bcftools:<br />
  file1=path to file1/${PROJECT_NAME}_Annotated.eff.vcf<br />
  file1=path to file2/${PROJECT_NAME}_Annotated.eff.vcf<br />
  bgzip $file1<br />
  bgzip $file1<br />
  file1=${file1}.gz<br />
  file2=${file2}.gz<br />
  tabix $file1<br />
  tabix $file2<br />
  bcftools merge -o {merged vcf name}.vcf $file1 $file2<br />
  grep -v "##" {merged vcf name}.vcf > {merged vcf name2}.vcf<br /> #remove header

    
## 3) Summarize FPKM, QC and cluster
   Input Materials:<br />
  -i {merged vcf name2}.vcf<br />
  
  Output Materials:<br />
  -o boxplot<br />
  -o hclust<br />
  -o Heatmap<br />

## 4). DESeq for gene level comparision
   Input Materials:<br />
  -i sample_table_zhunter.csv (table with samples annotation)<br />
  -i XX.bam.count (HTSeq-count files)<br />
  
  Output Materials:<br />
  -i gene expression analysis
  


 ## III. Requirement
  
  System requirements:
  >40 GB memory high-performance computing (HPC) environment, prefer Sun Grid Engine.
  Linux and Mac OS 64 bit system
  
  Software requirement:
  R 3.3.2<br />
  python2.7 in server
