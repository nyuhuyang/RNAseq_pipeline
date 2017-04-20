#!/bin/bash -l
#$ -N RNAmut
#$ -j y
#$ -m a
#$ -M yah2014@med.cornell.edu
#$ -l h_vmem=8G
#$ -pe smp 4
#$ -l zenodotus=TRUE
#$ -q *@@red



######################################################################
#---------------------Variables to be set-------------------------#
PROJECT_NAME="WTCHG"
path=/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/All_Projects_RNASeq_Data/ALL_Aligned_BAMS/BAMS_All_Samples
MutPath=/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_MUTS/${PROJECT_NAME}_MUTS
fa=/zenodotus/elementolab/scratch/akv3001/Indexed_genome/hg19_genome_index.fa
WM_bed=/home/yah2014/success_scripts/WM.bed
List_bam=/home/yah2014/success_scripts/SampleList_WTCHG_bam
SampleList_WTCHG=/home/yah2014/success_scripts/SampleList_WTCHG

SampleList="WTCHG_254364_705505 WTCHG_254364_705506 WTCHG_254364_706505 WTCHG_254364_706506 WTCHG_254364_707505 WTCHG_254364_707506 \
WTCHG_254364_708505 WTCHG_254364_708506 WTCHG_254364_709505 WTCHG_254364_709506 WTCHG_254364_710505 WTCHG_254364_710506 \
WTCHG_254364_711505 WTCHG_254364_711506 WTCHG_254364_712505 WTCHG_254365_705505 WTCHG_254365_705506 WTCHG_254365_706505 \
WTCHG_254365_706506 WTCHG_254365_707505 WTCHG_254365_707506 WTCHG_254365_708505 WTCHG_254365_709505 WTCHG_254365_710505 \
WTCHG_254365_711505 WTCHG_254365_712505 WTCHG_255088_705505 WTCHG_255088_705506 WTCHG_255088_706505 WTCHG_255088_706506 \
WTCHG_255088_707505 WTCHG_255088_707506 WTCHG_255088_708505 WTCHG_255088_708506 WTCHG_255088_709505 WTCHG_255088_709506 \
WTCHG_255088_710505 WTCHG_255088_710506 WTCHG_255088_711505 WTCHG_255088_711506 WTCHG_255088_712505 WTCHG_255089_705505 \
WTCHG_255089_705506 WTCHG_255089_706505 WTCHG_255089_706506 WTCHG_255089_707505 WTCHG_255089_707506 WTCHG_255089_708505 \
WTCHG_255089_709505 WTCHG_255089_710505 WTCHG_255089_711505 WTCHG_255089_712505 WTCHG_332092_706503 WTCHG_332092_706504 \
WTCHG_332092_707503 WTCHG_332092_707504 WTCHG_332092_708503 WTCHG_332092_708504 WTCHG_332092_709503 WTCHG_332092_709504 \
WTCHG_332092_710503 WTCHG_332092_710504 WTCHG_332092_711503 WTCHG_332092_711504 WTCHG_332092_712503 WTCHG_332092_712504 \
WTCHG_334784_706503 WTCHG_334784_706504 WTCHG_334784_707503 WTCHG_334784_707504 WTCHG_334784_708503 WTCHG_334784_708504 \
WTCHG_334784_709503 WTCHG_334784_709504 WTCHG_334784_710503 WTCHG_334784_710504 WTCHG_334784_711503 WTCHG_334784_711504 \
WTCHG_334784_712503 WTCHG_334784_712504"


echo "path="
echo $(ls -l $path)
echo "fa="
echo $(ls -l $fa)
echo " "
echo $(ls -l $MutPath)
echo " "
echo $(ls -l $WM_bed)
echo " "
echo $(ls -l $List_bam)
echo " "
######################################################################
#step 0  rsync---------------------------------------#
echo "####### rsync ##########"
cd $TMPDIR
for file in $SampleList; do
rsync -v -a --exclude 'Summary' $path/${file}_STAR/${file}_sorted.bam ./ #copy the sorted bam files
rsync -v -a --exclude 'Summary' $path/${file}_STAR/${file}_sorted.bam.bai ./ #copy the sorted bam.bai files
done
echo "BAMS PROJECT_NAMEs transfer Complished."
echo $(ls -l)

# step 1, Samtools mpileup--------------------------------
mkdir results
echo ""
echo "############# Samtools mpileup ##################"
CMD="/home/yah2014/Programs/samtools-1.4/samtools mpileup -f $fa -l $WM_bed -b $List_bam -o $TMPDIR/results/${PROJECT_NAME}.mpileup" 
#save mpileup PROJECT_NAME
eval $CMD

echo $(ls -l $TMPDIR/results/${PROJECT_NAME}.mpileup)

#---- step 2, Varscan--------------------------------
echo " "
echo "############# Varsan ##################"
/home/akv3001/jdk1.8.0_05/bin/java -jar $HOME/Programs/VarScan.v2.3.9.jar mpileup2cns $TMPDIR/results/${PROJECT_NAME}.mpileup \
--vcf-sample-list $SampleList_WTCHG --min-var-freq 0.01 --variants --output-vcf 1  > $TMPDIR/results/${PROJECT_NAME}_mutations.vcf
echo $(ls -l $TMPDIR/results/${PROJECT_NAME}_mutations.vcf)

#---- step 3, SnpEff annotation--------------------------------
echo ""
echo "############# SnpEff annotation ##################"
mkdir -p $TMPDIR/results/mutation_anno_verscan
cd $TMPDIR/results/mutation_anno_verscan
/home/akv3001/jdk1.8.0_05/bin/java -jar /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/snpEff/SnpSift.jar \
annotate /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/BED_GTF_References/dbsnp-V-146-All.vcf \
$TMPDIR/results/${PROJECT_NAME}_mutations.vcf  > $TMPDIR/results/mutation_anno_verscan/${PROJECT_NAME}_dbsnp_annotated.vcf
echo $(ls -l $TMPDIR/results/mutation_anno_verscan/${PROJECT_NAME}_dbsnp_annotated.vcf)

/home/akv3001/jdk1.8.0_05/bin/java -jar /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/snpEff/SnpSift.jar \
annotate /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/CosmicCodingMuts.vcf \
$TMPDIR/results/mutation_anno_verscan/${PROJECT_NAME}_dbsnp_annotated.vcf > $TMPDIR/results/mutation_anno_verscan/${PROJECT_NAME}_cosmic_dbsnp.vcf
echo $(ls -l $TMPDIR/results/mutation_anno_verscan/${PROJECT_NAME}_cosmic_dbsnp.vcf)

/home/akv3001/jdk1.8.0_05/bin/java -Xmx4G -jar  /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/snpEff/snpEff.jar \
-v  GRCh37.75 -c /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/snpEff/snpEff.config \
$TMPDIR/results/mutation_anno_verscan/${PROJECT_NAME}_cosmic_dbsnp.vcf  >$TMPDIR/results/mutation_anno_verscan/${PROJECT_NAME}_Annotated.eff.vcf
echo $(ls -l $TMPDIR/results/mutation_anno_verscan/)

mv snpEff_genes.txt ${PROJECT_NAME}_snpEff_genes.txt
mv snpEff_summary.html ${PROJECT_NAME}_snpEff_summary.html
echo $(ls -l $TMPDIR/results/mutation_anno_verscan/)
mkdir -p $TMPDIR/results/mutation_anno_verscan/${PROJECT_NAME}_results
mv * $TMPDIR/results/mutation_anno_verscan/${PROJECT_NAME}_results
#----------------Files Transfer---------------------------#
echo ""
echo "############# Files Transfer ##################"
rsync -rva $TMPDIR/results $MutPath
