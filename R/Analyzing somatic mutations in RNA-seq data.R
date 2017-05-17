#---
#   title: "Analyzing somatic mutations in RNA-seq data"
#   author: "Yang Hu"
#   date: "2017/5/3"
#   output:
#    html_document: default
#---
    

## Abstract

#It is part of a workflow for analyzing RNA-seq data from Waldenstrom Macroglobulinemia patients, in a cohort called "zhunter," from Harvard University. After Alignment, samtools sort, VarScan and Annotation, an annotated zhunter_Annotated.eff(data).vcf is generated like below:
    
#    |#CHROM	| POS	 | ID	    | REF	|ALT	|QUAL |FILTER |INFO       |FORMAT|ZH10_NWM07_CTTGTA_L005_R1_001              |...|
#    |:-----:|:------:|:--------:|:-----:|:-----:|:---:|:-----:|:---------:|:---:|:-----------------------------:|:-:|
#    |chr11	|47376915|.|G	    |C      |.	  |PASS   |..EFF=missense_variant..|GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR|0/0:22:11:11:11:0:0%:1E0:39:0:7:4:0:0||...|
#    
#    In part A, script will add gene name, label "MISSENSE/NONSENSE" , extract FREQ%, add counts, and filter out SNP like below:
#    
#    |CHROM	| POS	 |gene_POS	|gene	| ID| REF	|ALT	|QUAL |FILTER |INFO      |counts|ZH10_NWM07_CTTGTA_L005_R1_001|...|
#    |:-----:|:------:|:--------:|:-----:|:-----:|:---:|:-----:|:---------:|:---:|:-----------------------------:|:-:|:-----:|:------:|
#    |chr11	|47376915	|SPI1 47376915	|SPI1	|.	|G	|C	|.	|PASS	|MISSENSE	|9	|0|..|
#    
#    In part B, script will read FPKM results and cluster somatic mutations. 

## A) Reading Annotated vcf, add gene names and annotations
#### A-1) Setup enviroment
#detect OS and set enviroment
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/WTCHG/Mutation");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/WTCHG/Mutation");getwd();list.files()}
#Double check SampleList_zhunter_ID in sample_table_WTCHG.csv, must be the same order as WTCHG_Annotated.eff.vcf
sample_table_WTCHG<- read.csv("sample_table_WTCHG.csv",sep=",",header =T)
WTCHG_Sample_name <- as.data.frame(sample_table_WTCHG$Sample_name)
WTCHG_Sample_ID<- as.data.frame(sample_table_WTCHG$Sample_ID)
AnnotatedVcf_WTCHG <- read.csv("WTCHG_Annotated.eff.vcf",sep="\t",header =T) #Get WTCHG_Annotated.eff.vcf
dim(AnnotatedVcf_WTCHG)
#rename AnnotatedVcf_WTCHG
for (i in 1:nrow(WTCHG_Sample_name)){
    names(AnnotatedVcf_WTCHG)[9+i]<- as.character(WTCHG_Sample_name[i,1])#renmae after column 9
    }
#Generate a new file `AnnotatedVcf_zhunter_new.csv`.
AnnotatedVcf_zhunter_new <-AnnotatedVcf_zhunter[,1:2] 

#### A-2) Add gene names
AnnotatedVcf_zhunter_new[,"gene_POS"]<-NA #add a empty column to a dataframe 
colnames(AnnotatedVcf_zhunter_new)[1]<-"CHROM" # was #CHORM previously
for(i in 1:(dim(AnnotatedVcf_zhunter)[1]-1)){
    if (AnnotatedVcf_zhunter_new[i,"POS"] >= 27022522 && AnnotatedVcf_zhunter_new[i,"POS"] <=27108601 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr1")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("ARID1A ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "ARID1A"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 136871919 && AnnotatedVcf_zhunter_new[i,"POS"] <=136873813 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr2")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("CXCR4 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "CXCR4"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 38179969 && AnnotatedVcf_zhunter_new[i,"POS"] <=38184512 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr3")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("MYD88 ",AnnotatedVcf_zhunter_new[i,"POS"]);
    AnnotatedVcf_zhunter_new[i,"gene"] <- "MYD88"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 2743387 && AnnotatedVcf_zhunter_new[i,"POS"] <=2757752 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr4")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("TNIP2 ",AnnotatedVcf_zhunter_new[i,"POS"]);
    AnnotatedVcf_zhunter_new[i,"gene"] <- "TNIP2"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 78432907 && AnnotatedVcf_zhunter_new[i,"POS"] <=78532988 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr4")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("CXCL13 ",AnnotatedVcf_zhunter_new[i,"POS"]);
    AnnotatedVcf_zhunter_new[i,"gene"] <- "CXCL13"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 121613068 && AnnotatedVcf_zhunter_new[i,"POS"] <=121844021 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr4")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("PRDM5 ",AnnotatedVcf_zhunter_new[i,"POS"]);
    AnnotatedVcf_zhunter_new[i,"gene"] <- "PRDM5"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 122052564 && AnnotatedVcf_zhunter_new[i,"POS"] <=122137782 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr4")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("TNIP3 ",AnnotatedVcf_zhunter_new[i,"POS"]);
    AnnotatedVcf_zhunter_new[i,"gene"] <- "TNIP3"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 150409504 && AnnotatedVcf_zhunter_new[i,"POS"] <=150460645 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr5")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("TNIP1 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "TNIP1"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 143072604 && AnnotatedVcf_zhunter_new[i,"POS"] <=143266338 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr6")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("HIVEP2 ",AnnotatedVcf_zhunter_new[i,"POS"]);
    AnnotatedVcf_zhunter_new[i,"gene"] <- "HIVEP2"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 95947212 && AnnotatedVcf_zhunter_new[i,"POS"] <=96081655 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr22")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("IGLL5 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "IGLL5"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 38177969 && AnnotatedVcf_zhunter_new[i,"POS"] <=38186512 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr9")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("WNK2 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "WNK2"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 44953899 && AnnotatedVcf_zhunter_new[i,"POS"] <=44971759 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr11")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("TP53I11 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "TP53I11"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 47376409 && AnnotatedVcf_zhunter_new[i,"POS"] <=47400127 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr11")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("SPI1 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "SPI1"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 53773979 && AnnotatedVcf_zhunter_new[i,"POS"] <=53810226 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr12")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("SP1 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "SP1 "}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 43699412 && AnnotatedVcf_zhunter_new[i,"POS"] <=43785354 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr15")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("TP53BP1 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "TP53BP1"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 7571720 && AnnotatedVcf_zhunter_new[i,"POS"] <=7590868 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr17")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("TP53 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "TP53"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 20715727 && AnnotatedVcf_zhunter_new[i,"POS"] <=20840434 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr18")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("CABLES1 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "CABLES1"}
    else if(AnnotatedVcf_zhunter_new[i,"POS"] >= 23229960 && AnnotatedVcf_zhunter_new[i,"POS"] <=23238013 && AnnotatedVcf_zhunter_new[i,"CHROM"]=="chr22")
    {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0("IGLL5 ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- "IGLL5"}
    else {AnnotatedVcf_zhunter_new[i,"gene_POS"] <- paste0(i," ",AnnotatedVcf_zhunter_new[i,"POS"])
    AnnotatedVcf_zhunter_new[i,"gene"] <- NA}
}
AnnotatedVcf_zhunter<-AnnotatedVcf_zhunter[!is.na(AnnotatedVcf_zhunter_new[,"gene"]),]
AnnotatedVcf_zhunter_new<-AnnotatedVcf_zhunter_new[!is.na(AnnotatedVcf_zhunter_new[,"gene"]),]
rownames(AnnotatedVcf_zhunter_new)<-AnnotatedVcf_zhunter_new[,"gene_POS"]

#### A-3) Add MISSENSE & NONSSENSE annotations

AnnotatedVcf_zhunter_new<-cbind(AnnotatedVcf_zhunter_new,AnnotatedVcf_zhunter[,3:7])

AnnotatedVcf_zhunter_new[,"INFO"]<-NA #add a empty column to a dataframe 
AnnotatedVcf_zhunter_new[grep("MISSENSE",AnnotatedVcf_zhunter[,"INFO"]),"INFO"]<-"MISSENSE"  #Add MISSENSE INFO
AnnotatedVcf_zhunter_new[grep("NONSENSE",AnnotatedVcf_zhunter[,"INFO"]),"INFO"]<-"NONSENSE" #Add MISSENSE INFO
intersect <- intersect(grep("NONSENSE",AnnotatedVcf_zhunter[,"INFO"]),grep("MISSENSE",AnnotatedVcf_zhunter[,"INFO"])) #find intersect
AnnotatedVcf_zhunter_new[intersect,"INFO"]<-"MISSENSE & NONSENSE" #Add both


#### A-4) Extract FREQ% data from column 10 to 84


library(stringr)
AnnotatedVcf_zhunter_new[,"counts"]<-NA #add a empty column to a dataframe 
Sample<-AnnotatedVcf_zhunter[,c(10,11)] #Generate a temp dataframe with correct dimension
colnames(Sample)<-c("Sample1","Sample2")
for(i in 1:75){
    Sample[,"Sample1"]<-str_split_fixed(AnnotatedVcf_zhunter[,i+9], ":", 14)[,7] #Split column and only extract FREQ%
    Sample[,"Sample1"]<-as.numeric(sub("%","",Sample[,"Sample1"])) #convert character of percentage into numeric
    Sample[,"Sample1"][is.na(Sample[,"Sample1"])]<- 0 #replace NA values with zeros
    AnnotatedVcf_zhunter_new <- cbind(AnnotatedVcf_zhunter_new, Sample[,"Sample1"]) 
    colnames(AnnotatedVcf_zhunter_new)[i+11] <-as.character(SampleList_zhunter_ID[i,]) #rename columns
}

#### A-5) Add counts
#`Counts` are the total occurrence in 75 samples for one specific mutation.

counting <-function(x){#input is [1,85] vector
    c = 0 
    for(i in 1:(ncol(x)-11)){ #"Counts" is at column 11
        if(x[1,i+11]>0) # total occurrence of >0%
            c<-c+1
    }
    return(c)
}
for(i in 1:nrow(AnnotatedVcf_zhunter_new)){
    AnnotatedVcf_zhunter_new[i,"counts"]<-counting(AnnotatedVcf_zhunter_new[i,]) #
}

####  A-6) Filter out regular SNP
#    Generate file `AnnotatedVcf_zhunter_COSM.csv` contains COSM mutations only.
#    Generate file `AnnotatedVcf_zhunter_novel.csv` contains "COSM || MISSENSE || NONSENSE" mutations only.

AnnotatedVcf_zhunter_COSM <- AnnotatedVcf_zhunter_new[grep("COSM",AnnotatedVcf_zhunter_new[,"ID"]),]
AnnotatedVcf_zhunter_novel <- AnnotatedVcf_zhunter_new[AnnotatedVcf_zhunter_new[,"ID"] %in% ".",]
AnnotatedVcf_zhunter_novel<- AnnotatedVcf_zhunter_novel[c(grep("MISSENSE",AnnotatedVcf_zhunter_novel[,"INFO"]),
                                                          grep("NONSENSE",AnnotatedVcf_zhunter_novel[,"INFO"])),]
AnnotatedVcf_zhunter_novel<- rbind(AnnotatedVcf_zhunter_COSM,AnnotatedVcf_zhunter_novel)
AnnotatedVcf_zhunter_novel<- AnnotatedVcf_zhunter_novel[order(rownames(AnnotatedVcf_zhunter_novel)),]
AnnotatedVcf_zhunter_novel<- AnnotatedVcf_zhunter_novel[unique(AnnotatedVcf_zhunter_novel[,"gene_POS"]),]

### A-7) Run t-test (mutation type vs Healthy/Cancer Condition) to find out biomarker
sample_WTCHG.table <- read.csv("sample_table_WTCHG.csv")
Conditions<-sample_WTCHG.table$Condition

t_test <-data.frame()
for(i in 1:nrow(AnnotatedVcf_WTCHG_new)){
    t_test[i,1]<-rownames(AnnotatedVcf_WTCHG_new)[i]
    mut<-as.double(AnnotatedVcf_WTCHG_new[i,12:(11+nrow(WTCHG_Sample_ID))])
    t_test[i,2]<-t.test(mut ~ Conditions)[['p.value']]
    tmp<-t.test(mut ~ Conditions)[['estimate']]
    t_test[i,3]<-tmp[1]-tmp[2]
    t_test[i,4]<-AnnotatedVcf_WTCHG_new[i,"counts"] 
}
t_test<-t_test[order(t_test[,2]),]
colnames(t_test)<-c("Gene_POS","p.value","Cancer.mean - Healthy.mean","counts")
rownames(t_test)<-1:nrow(AnnotatedVcf_WTCHG_new)    

####  A-8) Export csv files


write.csv(AnnotatedVcf_zhunter_new,"AnnotatedVcf_zhunter_new.csv")
write.csv(AnnotatedVcf_zhunter_COSM,"AnnotatedVcf_zhunter_COSM.csv")
write.csv(AnnotatedVcf_zhunter_novel,"AnnotatedVcf_zhunter_novel.csv")
write.csv(t_test,"t_test.csv")

## B) Analyzing somatic mutations using FPKM
#### B-1) Get FPKM_zhunter Data
#       Read each `genes.FPKM_tracking` data into the file `FPKM_zhunter_temp.csv`.<br />
#       Generate file `FPKM_zhunter` contains FPKM data from 1st sample only.
#detect OS and set enviroment
if (Sys.info()[['sysname']]=="Darwin"){
        FPKM_WTCHG_files_path <-paste0("/Users/yah2014/Documents/Programs/R/FPKM/", 
                                       WTCHG_Sample_ID[,1],"_CuffLinks/genes.FPKM_tracking")
}
if(Sys.info()[['sysname']]=="Windows"){
        FPKM_WTCHG_files_path <-paste0("C:/Users/User/Documents/Programs/R/FPKM/", 
                                       WTCHG_Sample_ID[,1],"_CuffLinks/genes.FPKM_tracking")
}
    
FPKM_zhunter <-data.frame(x= str(0), y= integer(0)) #Generate empty FPKM_zhunter dataframe
FPKM_zhunter_temp  <- read.csv(FPKM_zhunter_files_path[1], header=T, sep="\t") #Read data from 1st sample
FPKM_zhunter_temp <- FPKM_zhunter_temp[order(FPKM_zhunter_temp$tracking_id),] #Reorder the gene name
FPKM_zhunter<-FPKM_zhunter_temp$FPKM #Fill up FPKM_zhunter dataframe with tracking id


#Generate file `FPKM_zhunter` contains FPKM data from all samples. This might take less than one minute.

for(i in 2:nrow(SampleList_zhunter)){ #skip the first one, which is added already
    FPKM_zhunter_temp  <- read.csv(FPKM_zhunter_files_path[i], header=T, sep="\t")
    FPKM_zhunter_temp <- FPKM_zhunter_temp[order(FPKM_zhunter_temp$tracking_id),] #Reorder the gene name
    FPKM_zhunter<-cbind(FPKM_zhunter,FPKM_zhunter_temp$FPKM)            
}


#Rename rows and columns of`FPKM_zhunter`

dim(FPKM_zhunter)
colnames(FPKM_zhunter) <-SampleList_zhunter_name$Sample_ID # or annotation$SimpleLabel
rownames(FPKM_zhunter) <-FPKM_zhunter_temp$tracking_id

 #set up cut off `FPKM_WTCHG`
FPKM_WTCHG<-FPKM_WTCHG[rowSums(FPKM_WTCHG)>1,] #set cut off >1

#### B-2) Charaterize somatic mutations


#Run quanlity control with boxplot

# set a margin
par(oma=c(3,3,3,3))  # all sides have 3 lines of space  
par(mar=c(2,2,2,2))
boxplot(log10(FPKM_WTCHG + 1),xlab="all test samples",cex=0.05, cex.axis=0.5, ylab = "log (base 10) RPKM + 1",main="Quanlity control for FPKM",las=2)
FPKM_WTCHG_clean <- FPKM_WTCHG[,!(colnames(FPKM_WTCHG) %in% c("ZH2_rnaWT2","ZH4_rnaWT4"))] #which sample to be excluded based on bad quanlity
boxplot(log10(FPKM_WTCHG_clean + 1),xlab="all test samples",cex=0.05, cex.axis=0.5, ylab = "log (base 10) RPKM + 1",main="Quanlity control for FPKM",las=2)

#Group sample with hclust

par(oma=c(3,3,3,3))
par(mar=c(2,2,2,2))
lm <-log(FPKM_zhunter_clean+1)
dm<-as.dist(1-cor(lm))
hm<-hclust(dm,method = "average")
p<-plot(hm,col = "black",cex = 0.5,cex.lab=0.75,main="Cluster Dendrogram for FPKM")


#Generate heatmap
# Prepare palette

my_palette1 <- colorRampPalette(c("antiquewhite2", "green", "blue"))(n = 100)
my_palette2 <- colorRampPalette(c("red", "green"))(n = 2)
my_palette3 <- colorRampPalette(c("yellow", "orange", "red","green","blue","purple"))(n = 6)


#rename

mut_gene <- AnnotatedVcf_WTCHG_novel[,12:(nrow(WTCHG_Sample_name)+11)]
for (i in 1:nrow(WTCHG_Sample_name)){
    names(mut_gene)[i]<- as.character(WTCHG_Sample_name[i,1])#renmae after column 9
}
#Clean up data
#mut_gene_clean <- mut_gene[,!(colnames(mut_gene) %in% c("ZH2_rnaWT2","ZH4_rnaWT4"))]



#Generate heatmap
library(gplots)
library(ComplexHeatmap)
par(oma=c(3,3,3,3))
par(mar=c(2,2,2,2)+1)

Heatmap(mut_gene_clean, 
        column_title = "Cluster of RNA expression vs DNA mutation genes in WM",
        show_row_names = TRUE,
        col=my_palette1, 
        row_title_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 2),
        column_names_gp = gpar(fontsize = 6),
        clustering_distance_columns = "euclidean",
        cluster_rows = FALSE,
        cluster_columns =  hm,
        row_dend_width = unit(5, "cm"),
        column_dend_height = unit(30, "mm"),
        heatmap_legend_param = list(title = "mutation rate %")
) 
