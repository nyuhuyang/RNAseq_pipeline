#list of server
aphrodite.med.cornell.edu
aristotle
pascal
paris
orleans
descartes
panda2.pbtech



#---------------------------#path of reference files------------------------------------------------
/home/akv3001/hg19_UCSC_ref.gtf         #akanksha's hg19 gtf
/home/akv3001/Mus_UCSC_ref.gtf          #akanksha's mice gtf
/zenodotus/elementolab/scratch/akv3001/Indexed_genome/hg19_Genomedir/      #akanksha's hg19 STARREF
/zenodotus/elementolab/scratch/akv3001/Indexed_genome/hg19_genome_index.fa


#---------------------------decompress-------------------------
#decompress bz2
bzip2 -d
#decompress tar
tar -xvf
#decompress gz
gunzip
#decompress tgz

#---------------------------#rsync command line------------------------------------------------

mac135103:Documents yah2014$ rsync -vra yah2014@aristotle.med.cornell.edu:/home/yah2014/success_scripts ./ #download script to Documents folder locally
mac135103:Documents yah2014$ rsync -va ~/Documents/*bam ssh -l aristotle.med.cornell.edu:/home/yah2014/ #upload script to home folder in aristotle server
mac135103:Desktop yah2014$ rsync -va Danwei ssh -l aphrodite.med.cornell.edu:/athena/elementolab/scratch/xfer/zhunter/All_Projects_RNASeq_Data #upload files to zhunter folder in aphrodite server

#---------------------------qsub-----
qlogin -l "h_vmem=40G, athena=true"
qstat -u '*'
qdel -j 
qsub -t 1-9 STAR-HTSeq_Danwei.sh
#---------------------------#path of storage folder------------------------------------------------

#athena zhunter folder: aphrodite-- yes
/athena/elementolab/scratch/xfer/zhunter/All_Projects_RNASeq_Data/
ALL_MUTS/
    WTCHG_MUTS/
        results/
    zhunter_MUTS/
zhunter/
#zenodotus folder1:
/zenodotus/elementolab/scratch/yah2014/

#zenodotus folder2:
/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/
ALL_Aligned_BAMS/
    BAMS_WTCHG/
    BAMS_zhunter/
ALL_COUNTS/
    WTCHG_Counts/
    zhunter_Counts/
ALL_CUFFLINKS/
    WTCHG_FPKM/
    zhunter_FPKM/
ALL_MUTS/
    WTCHG_MUTS/
        results/
    zhunter_MUTS/
        results/
/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_MUTS/zhunter_MUTS/results/


#---------------------------#How to use tools no at default path------------------------------------------------
/home/akv3001/STAR_2.4.0g1/bin/Linux_x86_64/STAR

/home/yah2014/Programs/samtools-1.4/samtools
/home/akv3001/Programs/samtools/samtools

/home/yah2014/.local/bin/python2.7
/home/yah2014/.local/bin/python2.7 /home/yah2014/.local/bin/htseq-count

/home/akv3001/jdk1.8.0_05/bin/java
/home/akv3001/jdk1.8.0_05/bin/java -Xmx4G -jar  /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mutation_tools/snpEff/snpEff.jar

/home/akv3001/jdk1.8.0_05/bin/java -jar $HOME/Programs/VarScan.v2.3.9.jar






#---------------------------#Install Python without root privilege------------------------------------------
ssh aristotle

mkdir -p /home/yah2014/.local/bin
cd /home/yah2014/.local/bin
wget --no-check-certificate https://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz
tar -xzf Python-2.7.11.tgz
cd Python-2.7.11
./configure --prefix=/home/yah2014/.local
make && make altinstall

#do a "which python" for the user and for root
#If you add this link, do a "which python" for the user and for root. If root is pointing to /usr/local/bin/python, remove the link you just added, and figure out something else.

#change PATH
PATH=/home/yah2014/.local/bin:$PATH
echo $PATH

#add PATH permanently
vi $HOME/.bashrc
export PATH=/home/yah2014/Programs/cellranger-2.0.0:$PATH
 

ln -s /home/yah2014/.local/bin/python2.7 /home/yah2014/.local/bin/python
# then type in python, it should be correct python I just installed.
ls -ltr /home/yah2014/.local/bin/python*
-rwxr-xr-x 1 yah2014 oelab 6271378 Mar  6 16:32 /home/yah2014/.local/bin/python2.7
-rwxr-xr-x 1 yah2014 oelab    1697 Mar  6 16:32 /home/yah2014/.local/bin/python2.7-config
lrwxrwxrwx 1 yah2014 oelab      34 Mar  6 16:40 /home/yah2014/.local/bin/python -> /home/yah2014/.local/bin/python2.7

#-----------------#How to install pip (python) to user without root access--------------------------------------
cd /home/yah2014/.local/bin
wget https://bootstrap.pypa.io/get-pip.py && python get-pip.py --user

pip install --user PySam


#-------------install htseq-count--------------------------------
pip install --user numpy
pip install --user htseq
pip install --user matplotlib
pip install --user hashlib
easy_install hashlib
easy_install micropython-hashlib
easy_install Tempfile
pip install --user redis_wrapped_mkstemp

#-------------install git--------------------------------
cd /home/yah2014/.local/bin
wget https://www.kernel.org/pub/software/scm/git/git-2.12.0.tar.gz
tar -xzf git-2.12.0.tar.gz
cd git-2.12.0
./configure --prefix=/home/yah2014/
make install

#---------------------------#Install R without root privilege------------------------------------------
# didn't work
wget --no-check-certificate https://cran.r-project.org/src/base/R-3/R-3.4.1.tar.gz
tar -xvfz R-3.4.1.tar.gz
cd R-3.4.1
./configure --prefix=/home/yah2014/Programs/R-3.4.1 --with-readline=no --with-x=no --enable-R-shlib LDFLAGS="-L/$HOME/Programs/zlib-1.2.11/lib -L/$HOME/Programs/bzip2-1.0.6/lib -L/$HOME/Programs/xz-5.2.2/lib -L/$HOME/Programs/pcre-8.40/lib -L/$HOME/Programs/curl-7.47.1/lib" CPPFLAGS="-I/$HOME/Programs/zlib-1.2.11/include -I/$HOME/Programs/bzip2-1.0.6/include -I/$HOME/Programs/xz-5.2.2/include -I/$HOME/Programs/pcre-8.40/include -I/$HOME/Programs/curl-7.47.1/include"
make
make check

#-----------------------------build R-devel without root privilege------------------------------------------
#http://pj.freefaculty.org/blog/?p=315
#especially to the conversation at buttom of page
cd ~/src
wget --no-check-certificate http://stat.ethz.ch/R/daily/R-devel_2017-06-29.tar.gz
tar xzvf R-devel_2017-06-29.tar.gz
cd ~/src/R-devel
mkdir builddir
cd ~/src/R-devel/builddir
../configure --prefix=$HOME/packages/R-devel --with-cairo --with-readline=no --with-x=no --with-jpeglib --with-tcltk --with-blas --enable-BLAS-shlib --with-lapack --enable-R-profiling '--enable-R-shlib' '--enable-memory-profiling'
#change --with-readline to --with-readline=no

make
/home/yah2014/src/R-devel/builddir/bin/R
export PATH=$HOME/src/R-devel/builddir/bin:$PATH
#---------------------------#Install zlib without root privilege-----
cd ~/src
wget --no-check-certificate http://www.zlib.net/zlib-1.2.11.tar.gz 
tar xvf zlib-1.2.11.tar.gz
cd zlib-1.2.11
./configure --prefix=$HOME/packages
make
make install
#4. Adjust the environment so R-devel builds will find packages installed there.
export PATH=$HOME/packages/bin:$PATH
export LD_LIBRARY_PATH=$HOME/packages/lib:$LD_LIBRARY_PATH 
export CFLAGS="-I$HOME/packages/include" 
export LDFLAGS="-L$HOME/packages/lib" 
 
#---------------------------#Install bzip2 without root privilege-----
cd ~/src
wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
tar xzvf bzip2-1.0.6.tar.gz
cd bzip2-1.0.6
#http://www.linuxfromscratch.org/lfs/view/development/chapter06/bzip2.html
vi Makefile
vi Makefile-libbz2_so
#changed the line CC =gcc to CC =gcc -fPIC

make -f Makefile-libbz2_so
 make clean
 make
 make -n install PREFIX=$HOME/packages
 make install PREFIX=$HOME/packages

#---------------------------#Install liblzma without root privilege-----
cd ~/src
wget http://tukaani.org/xz/xz-5.2.2.tar.gz
tar xzvf xz-5.2.2.tar.gz
cd xz-5.2.2
./configure --prefix=$HOME/packages
make -j3
make install

#---------------------------#Install pcre >= 8.10 without root privilege-----
cd ~/src
wget ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.40.tar.gz
tar xzvf pcre-8.40.tar.gz
cd pcre-8.40
./configure --enable-utf8 --prefix=$HOME/packages
make -j3
make install


#---------------------------#Install curl-7.53.1 without root privilege-----
cd ~/src
wget --no-check-certificate https://curl.haxx.se/download/curl-7.53.1.tar.gz
tar xzvf curl-7.53.1.tar.gz
cd curl-7.53.1
./configure --prefix=$HOME/packages
make -j3
make install

#---------------------------#Install R without root privilege------------------------------------------
cd ~/src/R-3.3.3
./configure --prefix=/home/yah2014/Programs/R-3.3.3 \
--with-readline=no --with-x=no --enable-R-shlib \
LDFLAGS="-L/$HOME/Programs/zlib-1.2.11/lib -L/$HOME/Programs/bzip2-1.0.6/lib \
-L/$HOME/Programs/xz-5.2.2/lib -L/$HOME/Programs/pcre-8.40/lib \
-L/$HOME/Programs/curl-7.47.1/lib" CPPFLAGS="-I/$HOME/Programs/zlib-1.2.11/include \
-I/$HOME/Programs/bzip2-1.0.6/include -I/$HOME/Programs/xz-5.2.2/include \
-I/$HOME/Programs/pcre-8.40/include -I/$HOME/Programs/curl-7.47.1/include"

#if fail, use
/home/akv3001/bin/R


#---------------------------#Install Rnaseqmut without root privilege------------------------------------------
cd ~/Programs
git clone https://github.com/davidliwei/rnaseqmut.git
cd rnaseqmut/bin
mv rnaseqmut.linux.x64 rnaseqmut
export PATH=$PATH:/home/yah2014/Programs/rnaseqmut/bin:$PATH

#---------------------------#export bin:$PATH------------------------------------------

export PATH=$PATH:/home/akv3001/Programs/R-3.2.1/bin:$PATH

#-------------------------------add header------------


/home/akv3001/Programs/samtools/samtools view -h ZH10_NWM07_CTTGTA_L005_R1_001_sorted.bam | head -10000 >ZH10_NWM07_CTTGTA_L005_R1_001_test.bam
/home/akv3001/Programs/samtools/samtools calmd -b ZH10_NWM07_CTTGTA_L005_R1_001_test.bam /zenodotus/elementolab/scratch/akv3001/Indexed_genome/hg19_genome_index.fa > ZH10_NWM07_CTTGTA_L005_R1_001_test.bam.MD

#---------------------------Install bcftools------------------------------------------
#wget https://github.com/samtools/bcftools/releases/download/1.4/bcftools-1.4.tar.bz2  didn't work in Linux
#try to make it in Mac  didn't work
#download to Mac, and rsync to Linux
bzip2 -d bcftools-1.4.tar.bz2
tar -xvf bcftools-1.4.tar
cd bcftools-1.4
make
make prefix=$HOME install

#---------------------------Install VarScan.v2.3.9.jar------------------------------------------

wget --no-check-certificate https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar/download?use_mirror=superb-dca2&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fvarscan%2F&use_mirror=superb-dca2




# step 0.5, Data preparation,recalculate MD tag
echo ""
echo "####### Step 0.5, Data preparation,recalculate MD tag  ##########"
for file in $BAMFILELIST; do
filebase=`basename $file`
CMD="samtools calmd -b $file /zenodotus/elementolab/scratch/akv3001/Indexed_genome/hg19_genome_index.fa > ${file}.MD"
echo "#### COMMAND LINE: $CMD"
eval $CMD
done
#---------grep XX from all files in the folder------------------------------------------------------------------
mv *1st.txt 1st
mv *rnaWT* rnaWT
path=/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_MUTS/zhunter_MUTS/results/1st/rnaWT
SampleList="echo $(ls $path)"
cd $path
for file in $SampleList; do
    if grep -q 38182641 "$file"; then
        echo $file
        grep 38182641 "$file"
    fi
done

############# filter mutation with gene name ########################################
cd mutation_anno

#------list interest mut gene names in file Mut_name  -------------------------
# Too many duplicate mutations at one position. Restart the process from zhunter_Annotated.eff2.vcf

grep MYD88 zhunter_Annotated.eff2.vcf | cut -f 1-7 > MYD88_name
:%s/chr3/chr3^IMYD88/g  #add one column named "MYD88"

grep IGLL5 zhunter_Annotated.eff2.vcf | cut -f 1-7 > IGLL5_name
:%s/chr22/chr22^IIGLL5/g

grep CXCR4 zhunter_Annotated.eff2.vcf | cut -f 1-7 > CXCR4_name
:%s/chr2/chr2^ICXCR4/g

grep ARID1A zhunter_Annotated.eff2.vcf | cut -f 1-7 > ARID1A_name
:%s/chr1/chr1^IARID1A/g

grep SPI1 zhunter_Annotated.eff2.vcf | cut -f 1-7 > SPI1_name
:%s/chr11/chr11^ISPI1/g
cat MYD88_name IGLL5_name CXCR4_name SPI1_name ARID1A_name | sort > Mut_counts.txt


head -1 ALLMUT.txt >head
vi head
:%s/^I/,/g
#copy head -1 ALLMUT.txt to Mut_name
:%s/,/^I/g

#------grep mut counts for each samples in file Mut_count  -------------------------

awk '{ print $3 > "pos_list"}' Mut_name.txt
cut -f 3 Mut_name.txt > pos_list
for pos in $(cat pos_list); do cat ALLMUT.txt|grep $pos;done >> Mut_count.txt

cut -f 1-4 Mut_count.txt
awk '{ ratio= ($11+$12)/($9+$10); print $1 $2 $3 $4 ratio > "pos_list"}' Mut_counts.test # ration =(reff+refv)/(altf+altv)

grep --line-buffered chr11 ALLMUT.txt |grep --line-buffered 47376544 |grep --line-buffered G | cat >> Mut_count.txt
for pos in $(cat pos_list); do cat zhunter_Annotated.eff2.vcf|grep $pos;done >> zhunter_Annotated.eff3.vcf
# Too many duplicate mutations at one position. Restart the process from Mut_counts.txt
vi Mut_counts.txt
:%s/;/^I/g7
:%s/,/^I/g
:%s/=/^I/g

cut -f 1-7 Mut_counts.txt >Mut_counts.name
# crate vcf contains chr pos rs ref alt ratio
awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t";for(i=8;i<=7+5*75;i=i+5) printf($i"\t"($(i+3)+$(i+4))/($(i+1)+$(i+2))"\t"); printf "\n"}' Mut_counts.txt> Mut_counts.counts

#awk: (FILENAME=Mut_counts.txt FNR=3) fatal: division by zero attempted
#test awk condition
awk '{for(i=8;i<=7+5*75;i=i+5) if($(i+1)+$(i+2)){printf($i"\t"($(i+3)+$(i+4))/($(i+1)+$(i+2))"\4/6/2017t")} else{printf 0"\t"}; printf "\n"}' Mut_counts.txt
awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t";for(i=8;i<=7+5*75;i=i+5) if($(i+1)+$(i+2)){printf($i"\t"($(i+3)+$(i+4))/($(i+1)+$(i+2))"\t")} else{printf 0"\t"}; printf "\n"}' Mut_counts.txt> Mut_counts.counts
#Some annotation has more than one column, groom manually
:%s/rs148489860^ICOSM3357306/COSM3357306/
:%s/rs387907272^ICOSM85940/COSM85940/
awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t";for(i=8;i<=7+5*75;i=i+5) if($(i+1)+$(i+2)){printf($i"\t"($(i+3)+$(i+4))/($(i+1)+$(i+2))"\t")} else{printf 0"\t"}; printf "\n"}' Mut_counts.txt> Mut_counts.counts
#check every sample name stays at the same position
awk '{for(i=8;i<=7+1*30;i=i+5) if($(i+1)+$(i+2)){printf($i"\t")} else{printf 0"\t"}; printf "\n"}' Mut_counts.txt |less
#Two many columns are not aligned. GIVE UP!

#debug
for (( c=1; c<=5; c++ ));do echo "Welcome $c times";done
for (( i=0;i<=5*75;i=i+5)); do echo "$i"; done
for i in 1 2 ; do awk '{ print $$i }' Mut_counts.test ;done # doesn't work
awk '{ ratio= ($11+$12)/($9+$10); printf $8 "\t" ratio "\n"}' Mut_counts.test # ration =(reff+refv)/(altf+altv)
awk '{for(i=8;i<=8+5*75;i=i+5) printf(($(i+3)+$(i+4))/($(i+1)+$(i+2))"\t")}' Mut_counts.txt >> Mut_counts.count #ration =(reff+refv)/(altf+altv)
paste Mut_counts.name Mut_counts.count > Mut_counts.counts
for (( c=1; c<=2; c++ ));do awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t";for(i=8;i<=7+5*75;i=i+5) printf($i"\t"($(i+3)+$(i+4))/($(i+1)+$(i+2))"\t"); printf "\n"}' Mut_counts.test |more;done   #print out first 7 columns and then sample name and mut ration =(reff+refv)/(altf+altv)

#-----------04/06/2017---------------

#count total number of files in directory linux
ls | wc -l  # ll | wc -l is wrong

# count total lines of all files in the folder
for file in $(ls); do wc -l $file; done

# count number of # in each vcf

for file in $(ls); do grep "#" $file |wc -l ; done   #=24

# remove # in each *_Annotated.eff.vcf > *_Annotated.eff.sim.vcf
for file in $(cat SampleList_zhunter); do grep -v "#" ${file}_results/${file}_Annotated.eff.vcf > ${file}_results/${file}_Annotated.eff.sim.vcf; done
# count number of mutation in each *_Annotated.eff.sim.vcf
for file in $(cat SampleList_zhunter); do wc -l ${file}_results/${file}_Annotated.eff.sim.vcf;done
# Summarize all mutations
for file in $(cat SampleList_zhunter); do cut -f 1-7 ${file}_results/${file}_Annotated.eff.sim.vcf |cat >> Mut_counts.name;done # remember use >> to create file

#combine all three lines, sort uniq -c > Mut_counts.name1 file contains all mutations happend in all samples including the number of occurrences
for file in $(cat SampleList_zhunter); do grep -v "#" ${file}_results/${file}_Annotated.eff.vcf | cut -f 1-5 |cat >> Mut_counts.name; done
sort Mut_counts.name |uniq -c >Mut_counts.name1

#only extract FREQ "Variant allele frequency" from column 8
# first separate all components in column 10 by replace : to tab,
# but first have to make sure don't replace : in other columns
#double check
grep ":" Mut_counts.name1
#remove annotation
cut --complement -s -f7-9 ZH10_NWM07_CTTGTA_L005_R1_001_results/ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.vcf| grep -v "#"> ZH10_NWM07_CTTGTA_L005_R1_001_results/ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf
for file in $(cat SampleList_zhunter); do grep -v "#" ${file}_results/${file}_Annotated.eff.vcf |cut --complement -s -f7-9 > ${file}_results/${file}_Annotated.eff.sim.vcf; done
#replace all : to tab
sed 's/:/\t/g' ZH10_NWM07_CTTGTA_L005_R1_001_results/ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf > ZH10_NWM07_CTTGTA_L005_R1_001_results/ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf1
cut -f 2,13 ZH10_NWM07_CTTGTA_L005_R1_001_results/ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf1 >ZH10_NWM07_CTTGTA_L005_R1_001_results/ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf2
#################combine all four steps##################------
for file in $(cat SampleList_zhunter); do grep -v "#" ${file}_results/${file}_Annotated.eff.vcf |cut --complement -s -f7-9| sed "s/:/\t/g" |cut -f 2,13 > ${file}_results/${file}_Annotated.eff.sim.vcf; done
----------###########################
# add file name to end column
#try sed
sed -i "s/$/\tZH10_NWM07_CTTGTA_L005_R1_001/" ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf
sed -i "s/$/\t$file"

for f in $(cat SampleList_zhunter); do sed -i "s/$/\t${f} " "${f}_results/${f}_Annotated.eff.sim.vcf"; done  # doesn't work, why???
for f in $(cat SampleList_zhunter); do sed -i "s/\$/\t${f}/g" ${f}_results/${f}_Annotated.eff.sim.vcf; done  #worked!!!!  use \$ and /g !!!!!!!
# try awk
for file in $(cat SampleList_zhunter); do awk '{printf $1"\t"$2"\t"${file}"\n"}' ${file}_results/${file}_Annotated.eff.sim.vcf > ${file}_results/${file}_Annotated.eff.sim.vcf1;done #work like shit....

# Swap two columns - awk
awk '{ t=$2; $2=$3; $3=t; print}' ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf
##########clear up files##############------
for file in $(cat SampleList_zhunter); do rm ${file}_results/${file}_Annotated.eff.sim.vcf; done
for file in $(cat SampleList_zhunter); do rm ${file}_results/${file}_Annotated.sim.vcf1; done
------############combine above two steps##################-----

for file in $(cat SampleList_zhunter); do sed "s/\$/\t${file}/g" ${file}_results/${file}_Annotated.eff.sim.vcf | awk '{ t=$2; $2=$3; $3=t; print$0}' |sed "s/ /\t/g">${file}_results/${file}_Annotated.eff.sim.vcf1; done
----###########################
#Asign all mutations from each vcf file to Mut_counts.name1
#first column is the key
# remove space generated by uniq -c
vi Mut_counts.name1
:%s/     //g
:%s/ /\t/g # change space to new line
# Swap two columns - awk
awk '{ t=$3; $3=$1; $1=t; print}' Mut_counts.name1 >Mut_counts.name2
#Asign all mutations from each vcf file to Mut_counts.name1
#--------
awk 'FNR==NR{a[$1]="\t"$2"\t"$3;next}{ print $0, a[$1]}' ZH10_NWM07_CTTGTA_L005_R1_001_results/ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf1 Mut_counts.name2 >Mut_counts.name3
#---------
#-----------04/07/2017---------------
# replce space into \t,remove mutiple \t\t, test it
for file in $(cat SampleList_zhunter); do sed "s/ /\t/g"  ${file}_Annotated.eff.sim.vcf2; done
for file in $(cat SampleList_zhunter); do sed "s/\t\t/\t/g"  ${file}_Annotated.eff.sim.vcf2; done
for file in $(cat SampleList_zhunter); do grep " " ${file}_Annotated.eff.sim.vcf2; done
# do all
###########################------
for file in $(cat SampleList_zhunter); do awk 'FNR==NR{a[$1]="\t"$2"\t"$3;next}{ print $0, a[$1]}' ${file}_results/${file}_Annotated.eff.sim.vcf1 Mut_counts.name2 | sed "s/ /\t/g" /dev/stdin |sed "s/\t\t/\t/g" /dev/stdin> ./eff.sim.vcf2/${file}_Annotated.eff.sim.vcf2;done
#----###########################
#test space and \t\t
for file in $(cat SampleList_zhunter); do grep " " ${file}_Annotated.eff.sim.vcf2; done

# align all vcf2
----###########################

# paste name with one vcf2
#remove -f1-7

cut -f8-9 ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf2 |paste Mut_counts.name3 /dev/stdin |less # works
cut -f8-9 ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf2 |paste Mut_counts.name3 /dev/stdin > Mut_counts.name3 #doesn't works
# above comman works for single vcf2, but will increse line number using multiple files.

paste Mut_counts.name3 <(cut -f8-9 ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf2) |less # works
paste Mut_counts.name3 <(cut -f8-9 ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf2) > Mut_counts.name3 #!doesn't works, only 2 column
paste Mut_counts.name3 <(cut -f8-9 ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf2) >> Mut_counts.name3 #! doesn't works, append

paste Mut_counts.name3 <(cut -f8-9 ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf2) |less
paste Mut_counts.name3 <(cut -f8-9 ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf2) > Mut_counts.name4
mv Mut_counts.name4 Mut_counts.name3
paste Mut_counts.name3 <(cut -f8-9 ZH10_rnaWT10_R1_.tmp_Annotated.eff.sim.vcf2) > Mut_counts.name4

# paste name with all vcf2 files
# To append text to a file you use >>. To overwrite the data currently in that file, you use >.
cp Mut_counts.name2 Mut_counts.name3

cut -f8-9 ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf2 > ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf3
paste Mut_counts.name3 ZH10_NWM07_CTTGTA_L005_R1_001_Annotated.eff.sim.vcf3 > Mut_counts.name4
mv Mut_counts.name4 Mut_counts.name3

cut -f8-9  ZH10_rnaWT10_R1_.tmp_Annotated.eff.sim.vcf2 > ZH10_rnaWT10_R1_.tmp_Annotated.eff.sim.vcf3
paste Mut_counts.name3 ZH10_rnaWT10_R1_.tmp_Annotated.eff.sim.vcf3 > Mut_counts.name4
mv Mut_counts.name4 Mut_counts.name3

#didn't work------------------------------------------------------------------------
for file in $(cat SampleList_zhunter)
do  cut -f8-9 ${file}_Annotated.eff.sim.vcf2 > ${file}_Annotated.eff.sim.vcf3
    paste Mut_counts.name3 ${file}_Annotated.eff.sim.vcf3 > Mut_counts.name3
done

for file in $(cat SampleList_zhunter)
do  cut -f8-9 ${file}_Annotated.eff.sim.vcf2 | paste Mut_counts.name3 /dev/stdin  > Mut_counts.name3 #have to switch the file name
done
#didn't work-------------------------------------------------------------

########################################################################-------
cp Mut_counts.name2 Mut_counts.name3
for file in $(cat SampleList_zhunter)
    do  cut -f8-9 ${file}_Annotated.eff.sim.vcf2 | paste Mut_counts.name3 /dev/stdin  > Mut_counts.name4
        mv Mut_counts.name4 Mut_counts.name3
done # It works. It didn't align correctly before since them have " " space in .sim.vcf2 files.
# add header
sed -i '1s/^/CHROM\tPOS\tcounts\tRS\tref\talt\t.\n/' Mut_counts.name2.txt
cut -f4-6 Mut_counts.name2.txt|less

for file in $(cat SampleList_zhunter);do sed -i '1s/^/CHROM\tPOS\tcounts\tRS\tref\talt\t.\t${file}\tFREQ\n/g' ${file}_Annotated.eff.sim.vcf2;done # add 1st line ${file} doesn't work '' doesn't work!!
for file in $(cat SampleList_zhunter);do sed -i "1s/^/CHROM\tPOS\tcounts\tRS\tref\talt\t.\t$file\t${file}_FREQ\n/g" ${file}_Annotated.eff.sim.vcf2;done # add 1st line
for file in $(cat SampleList_zhunter);do sed -i 1d ${file}_Annotated.eff.sim.vcf2;done # delete 1st line
for file in $(cat SampleList_zhunter);do mv ${file}_Annotated.eff.sim.vcf2 ${file}_Annotated.eff2.vcf;done
#Start R analysis

#-------########################

#2017/4/11 re-run Mpileup_Varscan_SnpEff_zhunter1------------------------------
#Since Akanksha suggest me to use merged mpileup file with all mutations from all files, but I realize I do mpileup seperately.
#So I decide to re-run Mpileup_Varscan_SnpEff_pipline.
#make --bam-list file contains List of input BAM files, one file per line
%s/ /_sorted.bam /g
%s/ /^M/g #replace space with new line ^M (control+V, control+M)
:%s/_R1_001/_R1_001_sorted.bam/g
:%s/_R1_.tmp/_R1_.tmp_sorted.bam/g


#2017/4/12---------------------
#make merged vcf file_zhunter from mpileup

#grep MISSENSE & COSMIC & .

#check each mutation with igv

#2017/4/13---------------------
#make heatmap
#Add phenotype in "INFO" column
mac135103:~ yah2014$ cd /Users/yah2014/Dropbox/Public/Olivier/R/zhunter/Expression

sync -rva yah2014@pascal.med.cornell.edu:/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_CUFFLINKS/zhunter_FPKM/ ./

#Start Mpileup_Varscan_SnpEff_WTCHG------
vi /home/yah2014/success_scripts/SampleList_zhunter_bam
:%s/\n/ /g #replace new line with sapce.

#2017/4/13---------------------
#make better heatmap
#replace long sample name
:%s/_L00.*_R1_001//g
############### work##################


############### work##################



#2017/4/16---------------------
#rerun Mpileup_Varscan_SnpEff_zhunter1.sh, add more genes
#remove "##" in zhunter_Annotated.eff20170416.vcf by vim or excel
#directly read zhunter_Annotated.eff20170416.vcf in to R
use http://atlasgeneticsoncology.org/ to find gene range

#2017/4/19---------------------
#merge two vcf files
install bcftools
install htslib
file1=/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_MUTS/zhunter_MUTS/results/mutation_anno_verscan_20170416/zhunter_results/zhunter_Annotated.eff20170416.vcf
file2=/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/yah2014/All_Projects_RNASeq_Data/ALL_MUTS/WTCHG_MUTS/results/mutation_anno_verscan/WTCHG_results/WTCHG_Annotated.eff.vcf
bgzip $file1
bgzip $file1
file1=${file1}.gz
file2=${file2}.gz
tabix $file1
tabix $file2

bcftools merge -o Annotated_WTCHG_zhunter_20170419.vcf $file1 $file2
grep -v "##" Annotated_WTCHG_zhunter_20170419.vcf > Annotated_WTCHG_zhunter_sim_20170419.vcf

#cluster fall apart
#will run PCA


#2017/5/4---------------------
#while run DESeqDataSetFromHTSeqCount
#if Error in `colnames<-`(`*tmp*`, value = 1:75) :
#attempt to set 'colnames' on an object with less than two dimensions
#clean the first row with the empty first column and "0" in the second column in HTSeq.count files
#use terminal, run below command lines to remove first line for all files
cd /Users/yah2014/Dropbox/Public/Olivier/R/ALL_COUNTS/zhunter_Counts
for file in $(ls);do echo "$(tail -n +2 $file)" > $file;done


R code is here:https://github.com/nyuhuyang/RNA-Seq-Analysis/blob/master/R/Analyzing%20somatic%20mutations%20in%20RNA-seq%20data.R
R Markdown demostration is here: https://github.com/nyuhuyang/RNA-Seq-Analysis/blob/master/vignettes/Analyzing%20somatic%20mutations%20in%20RNA-seq%20data.Rmd

#2017/05/16----------------------DESeq2 : three factor conditions--------------------------
#Deseq: Multiple Conditions Testing
If results is run without specifying contrast  or name , it will return the comparison of the last level of the last variable in the design formula over the first level of this variable. F
res_cancer <- results(dds,contrast = c("Condition","Cancer","Healthy"))

#! When Sample_name and Sample_ID are not highly replacable, Be very carefule when rename vcf file columns.
colnames(AnnotatedVcf_zhunter_new)[i+11] <-as.character(SampleList_zhunter_ID[i,]) #rename columns, here must be SampleList_zhunter_ID, not SampleList_zhunter_name.
#Double check SampleList_zhunter_ID in sample_table_WTCHG.csv, must be the same order as WTCHG_Annotated.eff.vcf

#2017/05/17---------------------t.test--------------
#What Biologist care most is raw data, then p-value and direction. Olivier said p-value doesn't mean anything if he can't see the trend from rawdata boxplot.
1. is p value, 2 is direction.
1.
t.test(mut ~ Conditions)$p.value
2.
t.test(mut ~ Conditions)$estimate[1]-t.test(mut ~ Conditions)$estimate[2]

Remember, Always make bigger/publishable Lengend and font.

#2017/05/18---------------------GSEA--------------
#Gene Annotation and pathway
#https://www.youtube.com/watch?v=vPRaUjiwwDQ
#Gene annotation website
#GeneCards (General information) www.genecards.org
#Uniport (protein level) www.uniprot.org
#Online Mendelian Inheritance in Man (gene history) www.omim.org
#NCBI Gene (NCBI resource) www.ncbi.nlm.nih.gov/gene

#Gene ontology

#Pathway database
#Kyoto Encyclopedia of Genes and Genomes www.genome.p/kegg/pathway.html
#Pathway interaction database pid.nci.nih.gov
#Reactome www.reactome.org , all on one

#Protein-protein interaction database
#ConsensusPathDB including compound

#David
#Gene_Ontology KEGG <---
#Problem for Functional Annotation Clustering
#1. Sampling issues
#2. Cut off bias
#3. Lost mild changes

#Gene set analysis
#https://www.youtube.com/watch?v=0qEqICa39U8
#Self-contained GSA  3/20 > 20/30000
#Competitive GSA -->GSEA

#2017/05/19---------------------GSEA--------------
source("http://bioconductor.org/biocLite.R")
## try http:// if https:// URLs are not supported

biocLite("BiocUpgrade")


#201/05/30----samtools--pileup--------
#didn't work if screen all variants. mpileup file become too large.


#2017/06/30
fastq --- head
sam, bam -----Samtools view

STAR if fastq is compressed
for gzipped files (*.gz) use --readFilesCommand zcat OR --readFilesCommand gunzip -c. 
for bzip2-compressed files, use --readFilesCommand bunzip2 -c

#Extract filename and extension in Bash
~% FILE="example.tar.gz"
~% echo "${FILE%%.*}"
example
~% echo "${FILE%.*}"
example.tar
~% echo "${FILE#*.}"
tar.gz
~% echo "${FILE##*.}"
gz

#$TMPDIR in qlogin
TMPDIR=/home/yah2014/TMPDIR
echo $(ls -l $TMPDIR/input/$gtf_name)

#check version
uname -a