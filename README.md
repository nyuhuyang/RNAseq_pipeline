----------------------------------------------<br />
RNA-Seq Analysis Pipeline: detecting mutations in RNA-Seq samples, exam genotype vs phenotype relationship<br />
Author: Yang Hu<br />
Department of Physiology and Biophysics<br />
Weill Cornell Medicine<br />
Email: yah2014@med.cornell.edu;Hu.Yang@nyu.edu<br />
----------------------------------------------<br />

## I. Introduction

This RNA-Seq Analysis pipeline has four sections:

1) <a href="https://github.com/nyuhuyang/RNAseq_pipeline/blob/master/ShellScript/Alignment_and_Counts.sh">Alignment and Counts</a></li>
Including STAR alignment, Samtools sort, HTSeq Count, and Cufflinks. Must be performed in the high-performance computing enviroment.

2) <a href="https://github.com/nyuhuyang/RNAseq_pipeline/blob/master/ShellScript/CallVar_and_Annotation.sh">Call Somatic mutations and annotate vcf</a></li>
Including Mpileup, Varscan, and SnpEff. Preferably be performed in the high-performance computing enviroment.

3) <a href="https://github.com/nyuhuyang/RNAseq_pipeline/blob/master/R/Analyzing_somatic_mutations.R">Summarize FPKM, QC and cluster</a></li>
Including mutation data clean, FPKM data clean, QC, hclust, and mutation/expression comparison. Can be performed on a local machine.

4) <a href="https://github.com/nyuhuyang/RNAseq_pipeline/blob/master/R/rnaseq_gene_level.R">DESeq for gene level comparision</a></li>
DESeq and downstram analysis. Can be performed on a local machine.


## II. Work flow

![plot of chunk Flow_work](vignettes/Flow_work.png)


### 1). Alignment_and_Counts.sh

  Input Materials:<br />
  -i `Sample.fastq.gz`<br />
  
  Output Materials:<br />
  -o samtools `Sample_name_sorted.bam`,`Sample_sorted.bam`,`Sample_sorted.bam.bai`<br />
  -o HTSeq `Sample.bam.count`<br />
  -o CuffLinks `genes.fpkm_tracking`,`genes.fpkm_tracking`, `genes.fpkm_tracking`, `genes.fpkm_tracking`<br />
  
  Linux bash shell script is here:<a href="https://github.com/nyuhuyang/RNAseq_pipeline/blob/master/ShellScript/Alignment_and_Counts.sh">Alignment and Counts</a></li>
  This bash shell is written for analyzing multiple projects. One project will have multiple RNA-seq samples.
  
  For example, there are total 155 Waldenström Macroglobulinemia (more than one terabyte) RNA-seq data from two projects("WTCHG" and "zhunter").<br />
  The PROJECT_NAME can be changed to fit different projects.
  
###  2). CallVar_and_Annotation.sh

  Input Materials:<br />
  -i `Sample_sorted.bam`,`Sample_sorted.bam.bai`<br />
  
  Output Materials:<br />
  -o samtools mpileup `${PROJECT_NAME}.mpileup`<br />
  -o VarScan `${PROJECT_NAME}_mutations.vcf`<br />
  -o snpEff `${PROJECT_NAME}_Annotated.eff.vcf`<br />
  
  Linux bash shell script is here:<a href="https://github.com/nyuhuyang/RNAseq_pipeline/blob/master/ShellScript/CallVar_and_Annotation.sh">Call Somatic mutations and annotate vcf</a></li>
  
### 2-1). If want to merge two Annotated.eff.vcf from two projects, run following script:


<!-- -->

    #merge multiple Annotated.eff.vcf files with bcftools:
    file1=path to file1/${PROJECT_NAME}_Annotated.eff.vcf
    file1=path to file2/${PROJECT_NAME}_Annotated.eff.vcf
    bgzip $file1
    bgzip $file1
    file1=${file1}.gz
    file2=${file2}.gz
    tabix $file1
    tabix $file2
    bcftools merge -o {merged vcf name}.vcf $file1 $file2
    grep -v "##" {merged vcf name}.vcf > {merged vcf name2}.vcf #remove header   
    
    
### 3) Summarize FPKM, QC and cluster
   Input Materials:<br />
  -i `{merged vcf name2}.vcf`<br />
  
  Output Materials:<br />
  -o boxplot<br />
  -o hclust<br />
  -o Heatmap<br />
  R code:<a href="https://github.com/nyuhuyang/RNAseq_pipeline/blob/master/R/Analyzing_somatic_mutations.R">Summarize FPKM, QC and cluster</a></li>
  R markdown: <a href="https://htmlpreview.github.io/?https://github.com/nyuhuyang/RNA-Seq-Analysis/blob/master/RMarkdown/Analyzing_somatic_mutations_in_RNA-seq_data.html">somatic variant analysis</a></li>"
  
### 4). DESeq for gene level comparision
   Input Materials:<br />
  -i `Sample_table_${PROJECT_NAME}.csv` (table with samples annotation)<br />
  -i `Sample.bam.count` (HTSeq-count files)<br />
  
  Output Materials:<br />
  -i gene expression analysis
  R code:<a href="https://github.com/nyuhuyang/RNAseq_pipeline/blob/master/R/rnaseq_gene_level.R">DESeq for gene level comparision</a></li>
  R markdown: <a href="https://htmlpreview.github.io/?https://github.com/nyuhuyang/RNA-Seq-Analysis/blob/master/RMarkdown/rnaseq_gene_level.html">Diferential expression analysis</a></li>"

  This is inspired by http://genomicsclass.github.io/book/pages/rnaseq_gene_level.html
 ## III. Requirement
  
  System requirements:
  >32 GB memory high-performance computing (HPC) environment, prefer Sun Grid Engine.
  Linux and Mac OS 64 bit system
  
  Software requirement:
  R > 3.3<br />
  python 2.7 in server
