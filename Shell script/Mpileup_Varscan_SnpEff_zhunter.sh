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
PROJECT_NAME="zhunter"
path=/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_Aligned_BAMS/BAMS_${PROJECT_NAME}
MutPath=/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_MUTS/${PROJECT_NAME}_MUTS
fa=/zenodotus/elementolab/scratch/akv3001/Indexed_genome/hg19_genome_index.fa
WM_bed=/home/yah2014/success_scripts/WM.bed
List_bam=/home/yah2014/success_scripts/SampleList_zhunter_bam
SampleList_zhunter=/home/yah2014/success_scripts/SampleList_zhunter
SampleList="ZH10_NWM07_CTTGTA_L005_R1_001 ZH10_rnaWT10_R1_.tmp ZH11_NWM09_ACAGTG_L005_R1_001 ZH11_rnaWT11_R1_.tmp \
ZH12_NWM10_GTGAAA_L006_R1_001 ZH12_rnaWT12_R1_.tmp ZH13_NWM11_GCCAAT_L006_R1_001 ZH13_rnaWT13_R1_.tmp \
ZH14_NWM12_CTTGTA_L007_R1_001 ZH14_rnaWT14_R1_.tmp ZH15_NWM13_ACAGTG_L007_R1_001 ZH15_rnaWT15_R1_.tmp \
ZH16_NWM15_GTGAAA_L008_R1_001 ZH16_rnaWT16_R1_.tmp ZH17_NWM16_GCCAAT_L008_R1_001 ZH17_rnaWT17_R1_.tmp \
ZH18_NWM17_CTTGTA_L001_R1_001 ZH18_rnaWT18_R1_.tmp ZH19_NWM18_ACAGTG_L001_R1_001 ZH19_rnaWT19_R1_.tmp \
ZH1_BCWM2_BM19plus_06_06_14_GATCAG_L001_R1_001 ZH1_NWM32_GCCAAT_L001_R1_001 ZH1_rnaWT1_R1_.tmp \
ZH20_NWM19_GTGAAA_L002_R1_001 ZH20_rnaWT20_R1_.tmp ZH21_NWM20_GCCAAT_L002_R1_001 ZH21_rnaWT21_R1_.tmp \
ZH22_NWM21_CTTGTA_L003_R1_001 ZH22_rnaWT22_R1_.tmp ZH23_NWM22_ACAGTG_L003_R1_001 ZH23_rnaWT23_R1_.tmp \
ZH24_NWM23_GTGAAA_L004_R1_001 ZH25_NWM24_GCCAAT_L004_R1_001 ZH26_NWM25_CTTGTA_L005_R1_001 ZH27_NWM26_ACAGTG_L005_R1_001 \
ZH28_NWM27_GTGAAA_L006_R1_001 ZH29_NWM28_GCCAAT_L006_R1_001 ZH2_BCWM2_P15_06_06_14_TAGCTT_L001_R1_001 \
ZH2_NWM33_CTTGTA_L001_R1_001 ZH2_rnaWT2_R1_.tmp ZH30_NWM29_CTTGTA_L007_R1_001 ZH31_NWM30_ACAGTG_L007_R1_001 \
ZH32_NWM31_GTGAAA_L004_R1_001 ZH33_WM26_GCCAAT_L004_R1_001 ZH34_WM27_CTTGTA_L005_R1_001 ZH35_WM28_ACAGTG_L005_R1_001 \
ZH36_WM29_GTGAAA_L006_R1_001 ZH37_WM01_GCCAAT_L006_R1_001 ZH38_WM02_CTTGTA_L007_R1_001 ZH39_WM03_ACAGTG_L007_R1_001 \
ZH3_NWM34_ACAGTG_L002_R1_001 ZH3_rnaWT3_R1_.tmp ZH40_WM04_GTGAAA_L008_R1_001 ZH41_WM05_GCCAAT_L008_R1_001 \
ZH42_WM06_CTTGTA_L001_R1_001 ZH43_WM07_ACAGTG_L001_R1_001 ZH44_WM08_GTGAAA_L002_R1_001 ZH45_WM09_GCCAAT_L002_R1_001 \
ZH46_WM10_CTTGTA_L003_R1_001 ZH4_rnaWT4_R1_.tmp ZH50_WM14_CTTGTA_L004_R1_001 ZH51_WM15_ACAGTG_L005_R1_001 \
ZH53_WM17_GCCAAT_L006_R1_001 ZH55_WM19_ACAGTG_L007_R1_001 ZH57_WM21_GCCAAT_L008_R1_001 ZH59_WM23_ACAGTG_L001_R1_001 \
ZH5_rnaWT5_R1_.tmp ZH62_BCWM1_CTTGTA_L002_R1_001 ZH63_MWCL1_ACAGTG_L002_R1_001 ZH6_rnaWT6_R1_.tmp \
ZH7_NWM04_ACAGTG_L003_R1_001 ZH7_rnaWT7_R1_.tmp ZH8_rnaWT8_R1_.tmp ZH9_NWM06_GCCAAT_L004_R1_001 ZH9_rnaWT9_R1_.tmp"

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
echo $(ls -l $SampleList_zhunter)
echo ""
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
--vcf-sample-list $SampleList_zhunter --min-var-freq 0.01 --variants --output-vcf 1  > $TMPDIR/results/${PROJECT_NAME}_mutations.vcf
echo $(ls -l $TMPDIR/results/${PROJECT_NAME}_mutations.vcf)

#---- step 3, SnpEff annotation--------------------------------
echo ""
echo "############# SnpEff annotation ##################"
mkdir -p $TMPDIR/results/mutation_anno_verscan_20170416
cd $TMPDIR/results/mutation_anno_verscan_20170416
/home/akv3001/jdk1.8.0_05/bin/java -jar /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/snpEff/SnpSift.jar \
annotate /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/BED_GTF_References/dbsnp-V-146-All.vcf \
$TMPDIR/results/${PROJECT_NAME}_mutations.vcf  > $TMPDIR/results/mutation_anno_verscan_20170416/${PROJECT_NAME}_dbsnp_annotated.vcf
echo $(ls -l $TMPDIR/results/mutation_anno_verscan_20170416/${PROJECT_NAME}_dbsnp_annotated.vcf)

/home/akv3001/jdk1.8.0_05/bin/java -jar /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/snpEff/SnpSift.jar \
annotate /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/CosmicCodingMuts.vcf \
$TMPDIR/results/mutation_anno_verscan_20170416/${PROJECT_NAME}_dbsnp_annotated.vcf > $TMPDIR/results/mutation_anno_verscan_20170416/${PROJECT_NAME}_cosmic_dbsnp.vcf
echo $(ls -l $TMPDIR/results/mutation_anno_verscan_20170416/${PROJECT_NAME}_cosmic_dbsnp.vcf)

/home/akv3001/jdk1.8.0_05/bin/java -Xmx4G -jar  /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/snpEff/snpEff.jar \
-v  GRCh37.75 -c /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/snpEff/snpEff.config \
$TMPDIR/results/mutation_anno_verscan_20170416/${PROJECT_NAME}_cosmic_dbsnp.vcf  \
>$TMPDIR/results/mutation_anno_verscan_20170416/${PROJECT_NAME}_Annotated.eff.vcf
echo $(ls -l $TMPDIR/results/mutation_anno_verscan_20170416/)

mv snpEff_genes.txt ${PROJECT_NAME}_snpEff_genes.txt
mv snpEff_summary.html ${PROJECT_NAME}_snpEff_summary.html
echo $(ls -l $TMPDIR/results/mutation_anno_verscan_20170416/)
mkdir -p $TMPDIR/results/mutation_anno_verscan_20170416/${PROJECT_NAME}_results
mv * $TMPDIR/results/mutation_anno_verscan_20170416/${PROJECT_NAME}_results
#----------------Files Transfer---------------------------#
echo ""
echo "############# Files Transfer ##################"
rsync -rva $TMPDIR/results $MutPath