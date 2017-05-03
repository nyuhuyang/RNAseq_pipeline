#!/bin/bash -l
#$ -N RNA-WM
#$ -j y
#$ -m a
#$ -M yah2014@med.cornell.edu
#$ -l h_vmem=8G
#$ -pe smp 4
#$ -q *@@red


#---------------------Variables to be set-------------------------#
PROJECT_NAME="zhunter"
path=/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/${PROJECT_NAME}
STARREF=/zenodotus/elementolab/scratch/akv3001/Indexed_genome/hg19_Genomedir/
gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/hg19_UCSC_ref_redone.gtf"
gtf_name=$(basename $gtf)
file=$(ls ${path} | tail -n +${SGE_TASK_ID}| head -1) # Uses job array for each sample in the folder
Sample=$(basename "$file")
Sample=${Sample%.*}  # remove .fastq
echo "path="
echo "$path
echo "STARREF=""
echo $(ls -l $STARREF)
echo "GTF="
echo $(ls -l $gtf)
echo $gtf_name
echo " "
echo "Sample=$Sample"

mkdir /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_Aligned_BAMS/BAMS_${PROJECT_NAME}
mkdir /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_CUFFLINKS/${PROJECT_NAME}_FPKM
mkdir /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_COUNTS/${PROJECT_NAME}_Counts

#----------------Files Transfer---------------------------#

cd $TMPDIR
mkdir input
echo "Start to transfer FASTQ files."
rsync -v -a -z --exclude 'Summary' $path/$file ./input
echo "total size"
echo $(ls -l ${TMPDIR}/input/$file)
echo "FASTQ files transfer Complished."
echo " "
echo "Start to transfer GTF file."
rsync -v -a $gtf ./input
echo $(ls -l ${TMPDIR}/input/$gtf_name)
echo "GTF transfer done"
echo " "

echo "Processing Sample"
echo " "

echo "Aligning and quantifying against Human-Hg19"


#--------Formatting multiple lane Data for STAR-------------------#
F1=$(ls $TMPDIR/input/$file/*R1*)
F2=$(ls $TMPDIR/input/$file/*R2*)
F1=$(echo $F1|tr ' ' ',') #Turn space into ,
F2=$(echo $F2|tr ' ' ',') #Turn space into ,
echo "F1= " ${F1}
echo "F2= " ${F2}

#-----------STAR Alignment Command--------------------------------#
echo "Processing Sample"
echo " "

mkdir ${TMPDIR}/${Sample}_STAR
cd ${TMPDIR}/${Sample}_STAR
echo "-------------------------------- "
echo 'Processing' ${Sample}
echo "Aligning FASTQ"
STAR --genomeDir ${STARREF} --readFilesIn ${F1} ${F2} --outSAMstrandField intronMotif --outFileNamePrefix ${Sample}  --runThreadN 4 --readFilesCommand zcat
echo "STAR alignment Complished"
echo "STAR output files:"
echo $(ls -l )
echo " "
#----------------Sorting SAM File by name-------------------
echo "Samtools start"
samtools view -bS *Aligned.out.sam > ${Sample}_Aligned.out.bam
samtools sort -n ${Sample}_Aligned.out.bam ${Sample}_name_sorted #for HTSeq
samtools sort ${Sample}_Aligned.out.bam ${Sample}_sorted # for Cufflinks

#----------------Samtools index ------------------------------

samtools index $TMPDIR/${Sample}_STAR/${Sample}_sorted.bam
rsync -r -v --exclude="*.sam" $TMPDIR/${Sample}_STAR /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_Aligned_BAMS/BAMS_${PROJECT_NAME}

echo $(ls -l )
#--------------Running HTSeq Count-----------------------------
# a Running htseq-count
echo " "
echo "HTSeq Count Start"
/home/yah2014/.local/bin/python2.7  /home/yah2014/.local/bin/htseq-count -s no -f bam $TMPDIR/${Sample}_STAR/${Sample}_name_sorted.bam $TMPDIR/input/$gtf_name  > $TMPDIR/${Sample}.bam.count

rsync -r -v $TMPDIR/${Sample}.bam.count  /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_COUNTS/${PROJECT_NAME}_Counts
echo "HTSeq Count End"

#-------------------CuffLinks Command---------------------------
echo " "
echo "Cufflinks Start"
cufflinks -q --max-bundle-frags 100000000  -p 2  -G $TMPDIR/input/$gtf_name -o $TMPDIR/${Sample}_CuffLinks  $TMPDIR/${Sample}_STAR/${Sample}_sorted.bam
echo "Cufflinks End"

#---------------------------------------------------------------
rsync -r -v $TMPDIR/${Sample}_CuffLinks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_CUFFLINKS/${PROJECT_NAME}_FPKM