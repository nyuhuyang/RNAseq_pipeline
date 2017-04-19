----------------------------------------------<br />
RNA-Seq Analysis Pipeline: detecting mutations in RNA-Seq samples, exam genotype vs phenotype relationship<br />
Author: Yang Hu<br />
Department of Physiology and Biophysics<br />
Weill Cornell Medicine<br />
Email: Hu.Yang@nyu.edu<br />
----------------------------------------------<br />

###0. Introduction

RNA-Seq Analysis Pipeline has three sections:

1) STAR-HTSeq.sh                    Linux bash shell script including STAR alignment, Samtools sort, HTSeq Count, and Cufflinks
2) Mpileup_Varscan_SnpEff.sh        Linux bash shell script including Mpileup, Varscan, and SnpEff
3) Annotatedeffvcf-FREQ_summary.R   R code including mutation data clean, FPKM data clean, QC, hclust, and mutation/expression comparison.

Starting Materials:
155 Waldenström Macroglobulinemia RNA-seq fastq files (more than 1 terabyte)
