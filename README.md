----------------------------------------------<br />
RNA-Seq Analysis Pipeline: detecting mutations in RNA-Seq samples, exam genotype vs phenotype relationship<br />
Author: Yang Hu<br />
Department of Physiology and Biophysics<br />
Weill Cornell Medicine<br />
Email: Hu.Yang@nyu.edu<br />
----------------------------------------------<br />

0. Introduction

RNA-Seq Analysis Pipeline has three sections:

1) STAR-HTSeq.sh                    Linux bash shell script including STAR alignment, Samtools sort, HTSeq Count, and Cufflinks
2) Mpileup_Varscan_SnpEff.sh        Linux bash shell script including Mpileup, Varscan, and SnpEff
3) Annotatedeffvcf-FREQ_summary.R   R code including mutation data clean, FPKM data clean, QC, hclust, and mutation/expression comparison.

1. Outline

Starting Materials:
RNA-seq fastq files

This pipeline is writen for 155 Waldenström Macroglobulinemia (more than 1 terabyte) RNA-seq data.

System requirements:
>40 GB memory high-performance computing (HPC) environment, perfer Sun Grid Engine.
Linux and Mac OS 64 bit system

Software requirement:
R 3.3.2
python2.7 in server
