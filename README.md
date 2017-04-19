----------------------------------------------<br />
RNA-Seq Analysis Pipeline: detecting mutations in RNA-Seq samples, exam genotype vs phenotype relationship<br />
Author: Yang Hu<br />
Department of Physiology and Biophysics<br />
Weill Cornell Medicine<br />
Email: Hu.Yang@nyu.edu<br />
----------------------------------------------<br />

I. Introduction

RNA-Seq Analysis Pipeline has three sections:

1) STAR-HTSeq.sh<br />
Linux bash shell script including STAR alignment, Samtools sort, HTSeq Count, and Cufflinks.

2) Mpileup_Varscan_SnpEff_WTCHG.sh & Mpileup_Varscan_SnpEff_zhunter.sh<br />
Lnux bash shell script including Mpileup, Varscan, and SnpEff.

3) Annotatedeffvcf-FREQ_summary.R<br />
R code including mutation data clean, FPKM data clean, QC, hclust, and mutation/expression comparison.


<br />
<br />
II. Outline
<br />
Starting Materials:
  RNA-seq fastq files

System requirements:
  >40 GB memory high-performance computing (HPC) environment, prefer Sun Grid Engine.
  Linux and Mac OS 64 bit system

  Software requirement:
  R 3.3.2<br />
  python2.7 in server
<br />
<br />
<br />
III. STAR-HTSeq.sh<br />
<br />
  Starting Materials:<br />
  RNA-seq fastq files<br />
  
  End Materials:<br />
  Aligned, sort, and indexed Bam files<br />
  HTSeq Counts<br />
  FPKM<br />
  
  This bash shell is written for analyzing multiple projects.
  
  There are total 155 Waldenström Macroglobulinemia (more than one terabyte) RNA-seq data.
  Those Waldenström Macroglobulinemia RNA-seq data are from two separate projects.<br />
  Project zhunter: 75 RNA-seq data from Harvard.<br />
  Project WTCHG: 80 RNA-seq data from France.<br />
  
  The STAR-HTSeq.sh is project specific. The current script is for project zhunter only.
  The PROJECT_NAME can be changed to fit different projects.
  
IV. Mpileup_Varscan_SnpEff.sh<br />
<br />

  Starting Materials:<br />
  Aligned, sort, and indexed Bam files<br />
  
  End Materials:<br />
  mpileup file<br />
  mutations.vcf<br />
  dbsnp_annotated.vcf<br />
  cosmic_dbsnp.vcf<br />
  Annotated.eff.vcf<br />
  snpEff_genes.txt<br />
  snpEff_summary.html<br />
  <br />
  <br />
V. Annotatedeffvcf-FREQ_summary.R
<br />
  Starting Materials:<br />
  Annotated.eff.vcf<br />
  
  End Materials:<br />
  boxplot<br />
  hclust<br />
  Heatmap<br />
