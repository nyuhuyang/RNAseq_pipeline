#---
#title: "Analyzing somatic mutations in RNA-seq data"
#author: "Yang Hu"
#date: "2017/5/3"
#output:
#html_document: default
#---

## Abstract

#It is part of a workflow for analyzing RNA-seq data from Waldenstrom Macroglobulinemia patients, in a cohort called "WTCHG," from Harvard University. After Alignment, samtools sort, VarScan and Annotation, an annotated WTCHG_Annotated.eff(data).vcf is generated like below:

#|#CHROM	| POS	 | ID	    | REF	|ALT	|QUAL |FILTER |INFO       |FORMAT|ZH10_NWM07_CTTGTA_L005_R1_001              |...|
#|:-----:|:------:|:--------:|:-----:|:-----:|:---:|:-----:|:---------:|:---:|:-----------------------------:|:-:|
#|chr11	|47376915|.|G	    |C      |.	  |PASS   |..EFF=missense_variant..|GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR|0/0:22:11:11:11:0:0%:1E0:39:0:7:4:0:0||...|

#In part A, script will add gene name, label "MISSENSE/NONSENSE" , extract FREQ%, add counts, and filter out SNP like below:

#|CHROM	| POS	 |gene_POS	|gene	| ID| REF	|ALT	|QUAL |FILTER |INFO      |counts|ZH10_NWM07_CTTGTA_L005_R1_001|...|
#|:-----:|:------:|:--------:|:-----:|:-----:|:---:|:-----:|:---------:|:---:|:-----------------------------:|:-:|:-----:|:------:|
#|chr11	|47376915	|SPI1 47376915	|SPI1	|.	|G	|C	|.	|PASS	|MISSENSE	|9	|0|..|

#In part B, script will read FPKM results and cluster somatic mutations.

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
#Generate a new file `AnnotatedVcf_WTCHG_new.csv`.

AnnotatedVcf_WTCHG_new <-AnnotatedVcf_WTCHG[,1:2]


#### A-2) Add gene names

AnnotatedVcf_WTCHG_new[,"gene_POS"]<-NA #add a empty column to a dataframe
colnames(AnnotatedVcf_WTCHG_new)[1]<-"CHROM" # was #CHORM previously
for(i in 1:(nrow(AnnotatedVcf_WTCHG)-1)){
    if (AnnotatedVcf_WTCHG_new[i,"POS"] >= 27022522 && AnnotatedVcf_WTCHG_new[i,"POS"] <=27108601 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr1")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("ARID1A ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "ARID1A"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 136871919 && AnnotatedVcf_WTCHG_new[i,"POS"] <=136873813 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr2")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("CXCR4 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "CXCR4"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 38179969 && AnnotatedVcf_WTCHG_new[i,"POS"] <=38184512 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr3")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("MYD88 ",AnnotatedVcf_WTCHG_new[i,"POS"]);
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "MYD88"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 2743387 && AnnotatedVcf_WTCHG_new[i,"POS"] <=2757752 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr4")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("TNIP2 ",AnnotatedVcf_WTCHG_new[i,"POS"]);
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "TNIP2"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 78432907 && AnnotatedVcf_WTCHG_new[i,"POS"] <=78532988 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr4")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("CXCL13 ",AnnotatedVcf_WTCHG_new[i,"POS"]);
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "CXCL13"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 121613068 && AnnotatedVcf_WTCHG_new[i,"POS"] <=121844021 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr4")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("PRDM5 ",AnnotatedVcf_WTCHG_new[i,"POS"]);
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "PRDM5"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 122052564 && AnnotatedVcf_WTCHG_new[i,"POS"] <=122137782 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr4")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("TNIP3 ",AnnotatedVcf_WTCHG_new[i,"POS"]);
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "TNIP3"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 150409504 && AnnotatedVcf_WTCHG_new[i,"POS"] <=150460645 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr5")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("TNIP1 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "TNIP1"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 143072604 && AnnotatedVcf_WTCHG_new[i,"POS"] <=143266338 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr6")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("HIVEP2 ",AnnotatedVcf_WTCHG_new[i,"POS"]);
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "HIVEP2"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 95947212 && AnnotatedVcf_WTCHG_new[i,"POS"] <=96081655 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr22")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("IGLL5 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "IGLL5"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 38177969 && AnnotatedVcf_WTCHG_new[i,"POS"] <=38186512 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr9")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("WNK2 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "WNK2"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 44953899 && AnnotatedVcf_WTCHG_new[i,"POS"] <=44971759 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr11")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("TP53I11 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "TP53I11"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 47376409 && AnnotatedVcf_WTCHG_new[i,"POS"] <=47400127 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr11")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("SPI1 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "SPI1"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 53773979 && AnnotatedVcf_WTCHG_new[i,"POS"] <=53810226 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr12")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("SP1 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "SP1 "}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 43699412 && AnnotatedVcf_WTCHG_new[i,"POS"] <=43785354 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr15")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("TP53BP1 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "TP53BP1"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 7571720 && AnnotatedVcf_WTCHG_new[i,"POS"] <=7590868 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr17")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("TP53 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "TP53"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 20715727 && AnnotatedVcf_WTCHG_new[i,"POS"] <=20840434 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr18")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("CABLES1 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "CABLES1"}
    else if(AnnotatedVcf_WTCHG_new[i,"POS"] >= 23229960 && AnnotatedVcf_WTCHG_new[i,"POS"] <=23238013 && AnnotatedVcf_WTCHG_new[i,"CHROM"]=="chr22")
    {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0("IGLL5 ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- "IGLL5"}
    else {AnnotatedVcf_WTCHG_new[i,"gene_POS"] <- paste0(" ",AnnotatedVcf_WTCHG_new[i,"POS"])
    AnnotatedVcf_WTCHG_new[i,"gene"] <- NA}
}
AnnotatedVcf_WTCHG<-AnnotatedVcf_WTCHG[!is.na(AnnotatedVcf_WTCHG_new[,"gene"]),]
AnnotatedVcf_WTCHG_new<-AnnotatedVcf_WTCHG_new[!is.na(AnnotatedVcf_WTCHG_new[,"gene"]),]
rownames(AnnotatedVcf_WTCHG_new)<-AnnotatedVcf_WTCHG_new[,"gene_POS"]

#### A-3) Add MISSENSE & NONSSENSE annotations

AnnotatedVcf_WTCHG_new<-cbind(AnnotatedVcf_WTCHG_new,AnnotatedVcf_WTCHG[,3:7])

AnnotatedVcf_WTCHG_new[,"INFO"]<-NA #add a empty column to a dataframe
AnnotatedVcf_WTCHG_new[grep("MISSENSE",AnnotatedVcf_WTCHG[,"INFO"]),"INFO"]<-"MISSENSE"  #Add MISSENSE INFO
AnnotatedVcf_WTCHG_new[grep("NONSENSE",AnnotatedVcf_WTCHG[,"INFO"]),"INFO"]<-"NONSENSE" #Add MISSENSE INFO
intersect <- intersect(grep("NONSENSE",AnnotatedVcf_WTCHG[,"INFO"]),grep("MISSENSE",AnnotatedVcf_WTCHG[,"INFO"])) #find intersect
AnnotatedVcf_WTCHG_new[intersect,"INFO"]<-"MISSENSE & NONSENSE" #Add both


#### A-4) Extract FREQ% data from column 10 to 84


library(stringr)
AnnotatedVcf_WTCHG_new[,"counts"]<-NA #add a empty column to a dataframe
Sample<-AnnotatedVcf_WTCHG[,c(10,11)] #Generate a temp dataframe with correct dimension
colnames(Sample)<-c("Sample1","Sample2")
for(i in 1:nrow(WTCHG_Sample_ID)){
    Sample[,"Sample1"]<-str_split_fixed(AnnotatedVcf_WTCHG[,i+9], ":", 14)[,7] #Split column and only extract FREQ%
    Sample[,"Sample1"]<-as.numeric(sub("%","",Sample[,"Sample1"])) #convert character of percentage into numeric
    Sample[,"Sample1"][is.na(Sample[,"Sample1"])]<- 0 #replace NA values with zeros
    AnnotatedVcf_WTCHG_new <- cbind(AnnotatedVcf_WTCHG_new, Sample[,"Sample1"])
    colnames(AnnotatedVcf_WTCHG_new)[i+11] <-as.character(WTCHG_Sample_ID[i,]) #rename columns
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
for(i in 1:nrow(AnnotatedVcf_WTCHG_new)){
    AnnotatedVcf_WTCHG_new[i,"counts"]<-counting(AnnotatedVcf_WTCHG_new[i,]) #
}

####  A-6) Filter out regular SNP
#Generate file `AnnotatedVcf_WTCHG_COSM.csv` contains COSM mutations only.<br />
#Generate file `AnnotatedVcf_WTCHG_novel.csv` contains "COSM || MISSENSE || NONSENSE" mutations only.

AnnotatedVcf_WTCHG_COSM <- AnnotatedVcf_WTCHG_new[grep("COSM",AnnotatedVcf_WTCHG_new[,"ID"]),]
AnnotatedVcf_WTCHG_novel <- AnnotatedVcf_WTCHG_new[AnnotatedVcf_WTCHG_new[,"ID"] %in% ".",]
AnnotatedVcf_WTCHG_novel<- AnnotatedVcf_WTCHG_novel[c(grep("MISSENSE",AnnotatedVcf_WTCHG_novel[,"INFO"]),
                                                          grep("NONSENSE",AnnotatedVcf_WTCHG_novel[,"INFO"])),]
AnnotatedVcf_WTCHG_novel<- rbind(AnnotatedVcf_WTCHG_COSM,AnnotatedVcf_WTCHG_novel)
AnnotatedVcf_WTCHG_novel<- AnnotatedVcf_WTCHG_novel[order(rownames(AnnotatedVcf_WTCHG_novel)),]
AnnotatedVcf_WTCHG_novel<- AnnotatedVcf_WTCHG_novel[unique(AnnotatedVcf_WTCHG_novel[,"gene_POS"]),]

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


write.csv(AnnotatedVcf_WTCHG_new,"AnnotatedVcf_WTCHG_new.csv")
write.csv(AnnotatedVcf_WTCHG_COSM,"AnnotatedVcf_WTCHG_COSM.csv")
write.csv(AnnotatedVcf_WTCHG_novel,"AnnotatedVcf_WTCHG_novel.csv")
write.csv(t_test,"t_test.csv")


## B) Analyzing somatic mutations using FPKM
#### B-1) Get FPKM_WTCHG Data
#Read each `genes.FPKM_tracking` data into the file `FPKM_WTCHG_temp.csv`.<br />
#Generate file `FPKM_WTCHG` contains FPKM data from 1st sample only.


#detect OS and set enviroment
if (Sys.info()[['sysname']]=="Darwin"){
        FPKM_WTCHG_files_path <-paste0("/Users/yah2014/Documents/Programs/R/FPKM/", 
                                       WTCHG_Sample_ID[,1],"_CuffLinks/genes.FPKM_tracking")
}
if(Sys.info()[['sysname']]=="Windows"){
        FPKM_WTCHG_files_path <-paste0("C:/Users/User/Documents/Programs/R/FPKM/", 
                                       WTCHG_Sample_ID[,1],"_CuffLinks/genes.FPKM_tracking")
}
FPKM_WTCHG <-data.frame(x= str(0), y= integer(0)) #Generate empty FPKM_WTCHG dataframe
FPKM_WTCHG_temp  <- read.csv(FPKM_WTCHG_files_path[1], header=T, sep="\t") #Read data from 1st sample
FPKM_WTCHG_temp <- FPKM_WTCHG_temp[order(FPKM_WTCHG_temp$tracking_id),] #Reorder the gene name
FPKM_WTCHG<-FPKM_WTCHG_temp$FPKM #Fill up FPKM_WTCHG dataframe with tracking id


#Generate file `FPKM_WTCHG` contains FPKM data from all samples. This might take less than one minute.

for(i in 2:nrow(WTCHG_Sample_name)){ #skip the first one, which is added already
    FPKM_WTCHG_temp  <- read.csv(FPKM_WTCHG_files_path[i], header=T, sep="\t")
    FPKM_WTCHG_temp <- FPKM_WTCHG_temp[order(FPKM_WTCHG_temp$tracking_id),] #Reorder the gene name
    FPKM_WTCHG<-cbind(FPKM_WTCHG,FPKM_WTCHG_temp$FPKM)
}


#Rename rows and columns of`FPKM_WTCHG`

dim(FPKM_WTCHG)
colnames(FPKM_WTCHG) <-WTCHG_Sample_name[,1] # or annotation$SimpleLabel
rownames(FPKM_WTCHG) <-FPKM_WTCHG_temp$tracking_id

        
#set up cut off `FPKM_WTCHG`
FPKM_WTCHG<-FPKM_WTCHG[rowSums(FPKM_WTCHG)>1,] #set cut off >1
dim(FPKM_WTCHG)
        
####  B-2) creat GCT file for GSEA
#creat "na" description
na_description<-data.frame(rep(NA,nrow(FPKM_WTCHG)))
suppressWarnings(FPKM_WTCHG.gct<-cbind(rownames(FPKM_WTCHG),na_description,FPKM_WTCHG))#Duplicated rownames will be lost during cbind
colnames(FPKM_WTCHG.gct)[1:2]<-c("NAME","Description")
write.csv(FPKM_WTCHG.gct,"FPKM_WTCHG.csv") # open file with excel, delete first column, save as tab delimited text file
####  B-3) creat cls file for GSEA
#need use terminal to %s/,/ /g

#### B-4) Charaterize somatic mutations

#Run quanlity control with boxplot
par(oma=c(3,3,3,3))  # all sides have 3 lines of space
par(mar=c(2,2,2,2))
boxplot(log10(FPKM_WTCHG + 1),xlab="all test samples",cex=0.05, cex.axis=0.5, ylab = "log (base 10) RPKM + 1",main="Quanlity control for FPKM",las=2)
FPKM_WTCHG_clean <- FPKM_WTCHG[,!(colnames(FPKM_WTCHG) %in% c("ZH2_rnaWT2","ZH4_rnaWT4"))] #which sample to be excluded based on bad quanlity
boxplot(log10(FPKM_WTCHG_clean + 1),xlab="all test samples",cex=0.05, cex.axis=0.5, ylab = "log (base 10) RPKM + 1",main="Quanlity control for FPKM",las=2)


#Group sample with hclust
par(oma=c(3,3,3,3))
par(mar=c(2,2,2,2))
lm <-log(FPKM_WTCHG_clean+1)
dm<-as.dist(1-cor(lm))
hm<-hclust(dm,method = "average")
p<-plot(hm,col = "black",cex = 0.75,main="Cluster Dendrogram for FPKM")


#Generate heatmap<br />
#Prepare palette

my_palette1 <- colorRampPalette(c("antiquewhite2", "green", "blue"))(n = 100)
my_palette2 <- colorRampPalette(c("red", "green"))(n = 2)
my_palette3 <- colorRampPalette(c("yellow", "orange", "red","green","blue","purple"))(n = 6)



#rename

mut_gene <- AnnotatedVcf_WTCHG_novel[,12:(nrow(WTCHG_Sample_name)+11)]
for (i in 1:nrow(WTCHG_Sample_name)){
    names(mut_gene)[i]<- as.character(WTCHG_Sample_name[i,1])#renmae after column 9
}
#mut_gene_clean <- mut_gene[,!(colnames(mut_gene) %in% c("ZH2_rnaWT2","ZH4_rnaWT4"))]



#Generate heatmap

library(gplots)
library(ComplexHeatmap)
par(oma=c(3,3,3,3))
par(mar=c(2,2,2,2)+1)

Heatmap(mut_gene,
        column_title = "Cluster of RNA expression vs DNA mutation genes in WM",
        show_row_names = TRUE,
        col=my_palette1,
        row_title_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 8),
        clustering_distance_columns = "euclidean",
        cluster_rows = FALSE,
        cluster_columns =  hm,
        row_dend_width = unit(5, "cm"),
        column_dend_height = unit(30, "mm"),
        heatmap_legend_param = list(title = "mutation rate %")
)


+
Heatmap(Condition,
        show_row_names = TRUE,
        col=my_palette2, 
        width = unit(1, "cm"),
        heatmap_legend_param = list(title = "Disease")
) +
    Heatmap(Population,
            show_row_names = FALSE,
            col=my_palette3, 
            width = unit(1, "cm"),
            heatmap_legend_param = list(title = "Population")
    )
