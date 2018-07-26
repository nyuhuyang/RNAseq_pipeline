#!/bin/bash -l
#$ -N RNAmut
#$ -j y
#$ -m a
#$ -M yah2014@med.cornell.edu
#$ -l h_vmem=2G
#$ -l zenodotus=TRUE
#$ -q *@@red
######################################################################
#---------------------Variables to be set-------------------------#
path=/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_Aligned_BAMS/BAMS_zhunter
PROJECT_NAME="zhunter"
fa='/zenodotus/elementolab/scratch/akv3001/Indexed_genome/hg19_genome_index.fa'
echo "path="
echo "$path"
echo "fa="
echo "$fa"
echo " "
reads_file=$(tail -n +${SGE_TASK_ID} $reads_list| head -1)
base_file=${reads_file%_sorted.bam}
echo $(ls -l $path/${base_file}_STAR/$reads_file)

#mkdir /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_MUTS/${PROJECT_NAME}_MUTS
######################################################################
#------------rsync---------------------------------------#
cd $TMPDIR
echo "Start to transfer BAM file."
#for Sample in $SampleList; do
rsync -v -a --exclude 'Summary' $path/${base_file}_STAR/$reads_file ./ #copy the sorted bam files
echo "Start to transfer fa file."
rsync -v -a --exclude 'Summary' $fa ./ #copy the sorted bam files
echo "BAMS and FA files transfer Complished."
echo $(ls -l)

#####################################################################
#----------------Samtools recalculate MD tag and index ------------------------------
#for Sample in $SampleList; do
/home/yah2014/Programs/samtools-1.4/samtools calmd -b ${base_file}_sorted.bam hg19_genome_index.fa > ${base_file}_sorted.MD.bam
/home/yah2014/Programs/samtools-1.4/samtools index ${base_file}_sorted.MD.bam
echo $(ls -l)
####################
#----------------Files Transfer---------------------------#
echo "Start to transfer BAM file."
rsync -r -v ./*.MD.* /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_MUTS/${PROJECT_NAME}_MUTS
