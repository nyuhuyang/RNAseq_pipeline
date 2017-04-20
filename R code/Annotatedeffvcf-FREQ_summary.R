library(stringr)
library(gplots)
library(ComplexHeatmap)
# Clean Annotated.eff.vcf file
# Orignal file is like below:
#----------------------------------

# #CHROM	POS	ID	REF	ALT	Sample1    Sample2  Sample3
#[1,] chr1	27092944	.	A	G   0/1:32:105:105:9:9:50%:5.1645E-4:38:21:4:5:9:0	./.:.:47	1/1:57:86:84:3:13:81.25%:1.6121E-6:37:17:1:2:10:3

# What I want is like below:
# #CHROM	POS	ID	REF	ALT	ZH10_NWM07_CTTGTA_L005_R1_001    ZH10_rnaWT10_R1_.tmp   ZH11_rnaWT11_R1_.tmp
#[1,] chr1	27092944	.	A	G   50%        81.25%
#Still keep the INFO column.

#Step 0
#Remove ## header manually , or in linux grep -v "##" zhunter_Annotated.eff.vcf  
#Replace Sample names.

########################################################################
### A) Reading Annotated vcf and transform it into matrix format, add gene name and annotation
########################################################################

setwd("/Users/yah2014/Dropbox/Public/Olivier/R/WTCHG/Mutation");getwd();list.files() #set_up_enviroment
SampleList <- read.csv("SampleList_WTCHG", header=T)#Get_SampleList_Data-
SampleList_condition<- read.csv("SampleList_WTCHG_condition", header=T)
AnnotatedVcf <- read.csv("WTCHG_Annotated.eff.vcf",sep="\t",header =T) #Get zhunter_Annotated.eff.vcf-
dim(AnnotatedVcf)
AnnotatedVcf_new <-AnnotatedVcf[,1:2] #generate a new matrix to keep Chrom and pos

#-------Add gene names----------
AnnotatedVcf_new[,"gene_POS"]<-NA #add a empty column to a dataframe 
colnames(AnnotatedVcf_new)[1]<-"CHROM" # was #CHORM previously
for(i in 1:(nrow(AnnotatedVcf)-1)){
    if (AnnotatedVcf_new[i,"POS"] >= 27022522 && AnnotatedVcf_new[i,"POS"] <=27108601 && AnnotatedVcf_new[i,"CHROM"]=="chr1")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("ARID1A ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "ARID1A"}
    else if(AnnotatedVcf_new[i,"POS"] >= 136871919 && AnnotatedVcf_new[i,"POS"] <=136873813 && AnnotatedVcf_new[i,"CHROM"]=="chr2")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("CXCR4 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "CXCR4"}
    else if(AnnotatedVcf_new[i,"POS"] >= 38179969 && AnnotatedVcf_new[i,"POS"] <=38184512 && AnnotatedVcf_new[i,"CHROM"]=="chr3")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("MYD88 ",AnnotatedVcf_new[i,"POS"]);
        AnnotatedVcf_new[i,"gene"] <- "MYD88"}
    else if(AnnotatedVcf_new[i,"POS"] >= 2743387 && AnnotatedVcf_new[i,"POS"] <=2757752 && AnnotatedVcf_new[i,"CHROM"]=="chr4")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("TNIP2 ",AnnotatedVcf_new[i,"POS"]);
        AnnotatedVcf_new[i,"gene"] <- "TNIP2"}
    else if(AnnotatedVcf_new[i,"POS"] >= 78432907 && AnnotatedVcf_new[i,"POS"] <=78532988 && AnnotatedVcf_new[i,"CHROM"]=="chr4")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("CXCL13 ",AnnotatedVcf_new[i,"POS"]);
        AnnotatedVcf_new[i,"gene"] <- "CXCL13"}
    else if(AnnotatedVcf_new[i,"POS"] >= 121613068 && AnnotatedVcf_new[i,"POS"] <=121844021 && AnnotatedVcf_new[i,"CHROM"]=="chr4")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("PRDM5 ",AnnotatedVcf_new[i,"POS"]);
        AnnotatedVcf_new[i,"gene"] <- "PRDM5"}
    else if(AnnotatedVcf_new[i,"POS"] >= 122052564 && AnnotatedVcf_new[i,"POS"] <=122137782 && AnnotatedVcf_new[i,"CHROM"]=="chr4")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("TNIP3 ",AnnotatedVcf_new[i,"POS"]);
        AnnotatedVcf_new[i,"gene"] <- "TNIP3"}
    else if(AnnotatedVcf_new[i,"POS"] >= 150409504 && AnnotatedVcf_new[i,"POS"] <=150460645 && AnnotatedVcf_new[i,"CHROM"]=="chr5")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("TNIP1 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "TNIP1"}
    else if(AnnotatedVcf_new[i,"POS"] >= 143072604 && AnnotatedVcf_new[i,"POS"] <=143266338 && AnnotatedVcf_new[i,"CHROM"]=="chr6")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("HIVEP2 ",AnnotatedVcf_new[i,"POS"]);
        AnnotatedVcf_new[i,"gene"] <- "HIVEP2"}
    else if(AnnotatedVcf_new[i,"POS"] >= 95947212 && AnnotatedVcf_new[i,"POS"] <=96081655 && AnnotatedVcf_new[i,"CHROM"]=="chr22")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("IGLL5 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "IGLL5"}
    else if(AnnotatedVcf_new[i,"POS"] >= 38177969 && AnnotatedVcf_new[i,"POS"] <=38186512 && AnnotatedVcf_new[i,"CHROM"]=="chr9")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("WNK2 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "WNK2"}
    else if(AnnotatedVcf_new[i,"POS"] >= 44953899 && AnnotatedVcf_new[i,"POS"] <=44971759 && AnnotatedVcf_new[i,"CHROM"]=="chr11")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("TP53I11 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "TP53I11"}
    else if(AnnotatedVcf_new[i,"POS"] >= 47376409 && AnnotatedVcf_new[i,"POS"] <=47400127 && AnnotatedVcf_new[i,"CHROM"]=="chr11")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("SPI1 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "SPI1"}
    else if(AnnotatedVcf_new[i,"POS"] >= 53773979 && AnnotatedVcf_new[i,"POS"] <=53810226 && AnnotatedVcf_new[i,"CHROM"]=="chr12")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("SP1 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "SP1 "}
    else if(AnnotatedVcf_new[i,"POS"] >= 43699412 && AnnotatedVcf_new[i,"POS"] <=43785354 && AnnotatedVcf_new[i,"CHROM"]=="chr15")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("TP53BP1 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "TP53BP1"}
    else if(AnnotatedVcf_new[i,"POS"] >= 7571720 && AnnotatedVcf_new[i,"POS"] <=7590868 && AnnotatedVcf_new[i,"CHROM"]=="chr17")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("TP53 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "TP53"}
    else if(AnnotatedVcf_new[i,"POS"] >= 20715727 && AnnotatedVcf_new[i,"POS"] <=20840434 && AnnotatedVcf_new[i,"CHROM"]=="chr18")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("CABLES1 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "CABLES1"}
    else if(AnnotatedVcf_new[i,"POS"] >= 23229960 && AnnotatedVcf_new[i,"POS"] <=23238013 && AnnotatedVcf_new[i,"CHROM"]=="chr22")
        {AnnotatedVcf_new[i,"gene_POS"] <- paste0("IGLL5 ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- "IGLL5"}
    else {AnnotatedVcf_new[i,"gene_POS"] <- paste0(i," ",AnnotatedVcf_new[i,"POS"])
        AnnotatedVcf_new[i,"gene"] <- NA}
}
AnnotatedVcf<-AnnotatedVcf[!is.na(AnnotatedVcf_new[,"gene"]),]
AnnotatedVcf_new<-AnnotatedVcf_new[!is.na(AnnotatedVcf_new[,"gene"]),]

rownames(AnnotatedVcf_new)<-AnnotatedVcf_new[,"gene_POS"]
#-------Add MISSENSE & NONSSENSE INFO----------
AnnotatedVcf_new<-cbind(AnnotatedVcf_new,AnnotatedVcf[,3:7])

AnnotatedVcf_new[,"INFO"]<-NA #add a empty column to a dataframe 
AnnotatedVcf_new[grep("MISSENSE",AnnotatedVcf[,"INFO"]),"INFO"]<-"MISSENSE"  #Add MISSENSE INFO
AnnotatedVcf_new[grep("NONSENSE",AnnotatedVcf[,"INFO"]),"INFO"]<-"NONSENSE" #Add MISSENSE INFO
intersect <- intersect(grep("NONSENSE",AnnotatedVcf[,"INFO"]),grep("MISSENSE",AnnotatedVcf[,"INFO"])) #find intersect
AnnotatedVcf_new[intersect,"INFO"]<-"MISSENSE & NONSENSE" #Add both
AnnotatedVcf_new[,"counts"]<-NA #add a empty column to a dataframe 

##################################################################
### B) extract FREQ% only from rest columns, clean data
##################################################################

Sample<-AnnotatedVcf[,c(10,11)] #Generate a temp dataframe with correct dimension
colnames(Sample)<-c("Sample1","Sample2")
for(i in 1:nrow(SampleList_condition)){
    Sample[,"Sample1"]<-str_split_fixed(AnnotatedVcf[,i+9], ":", 14)[,7] #Split column and only extract FREQ%
    Sample[,"Sample1"]<-as.numeric(sub("%","",Sample[,"Sample1"])) #convert character of percentage into numeric
    Sample[,"Sample1"][is.na(Sample[,"Sample1"])]<- 0 #replace NA values with zeros
    AnnotatedVcf_new <- cbind(AnnotatedVcf_new, Sample[,"Sample1"]) 
    colnames(AnnotatedVcf_new)[i+11] <-as.character(SampleList_condition[i,]) #rename columns
}
#-------Add counts, what's the total occurrence in 75 samples for one specific mutation----------
counting <-function(x){#input is [1,85] vector
    c = 0 
    for(i in 12:86){
        if(x[1,i]>0) # total occurrence of >0%
            c<-c+1
        }
    return(c)
    }
for(i in 1:(dim(AnnotatedVcf_new)[1]-1)){
    AnnotatedVcf_new[i,"counts"]<-counting(AnnotatedVcf_new[i,]) #
}

#--------grep COSM || MISSENSE || NONSENSE only mutations--------------
AnnotatedVcf_COSM <- AnnotatedVcf_new[grep("COSM",AnnotatedVcf_new[,"ID"]),]
AnnotatedVcf_novel <- AnnotatedVcf_new[AnnotatedVcf_new[,"ID"] %in% ".",]
AnnotatedVcf_novel<- AnnotatedVcf_novel[c(grep("MISSENSE",AnnotatedVcf_novel[,"INFO"]),
                                        grep("NONSENSE",AnnotatedVcf_novel[,"INFO"])),]
AnnotatedVcf_novel<- rbind(AnnotatedVcf_COSM,AnnotatedVcf_novel)
AnnotatedVcf_novel<- AnnotatedVcf_novel[order(rownames(AnnotatedVcf_novel)),]
AnnotatedVcf_novel<- AnnotatedVcf_novel[unique(AnnotatedVcf_novel[,"gene_POS"]),]

#--------Add phenotype in "INFO" column--------------
#/// AnnotatedVcf_COSM <- AnnotatedVcf_COSM[-c(5, 9, 15, 16), ]
#/// AnnotatedVcf_COSM[,"INFO"]<-c("MISSENSE; P->L",
#///                              "NONSENSE; Q->STOP",
#///                              "NA; frameshift",
#///                              "MISSENSE: R->C",
#///                              "MISSENSE; M->K",
#///                              "MISSENSE; V->G",
#///                              "MISSENSE; A->T",
#///                              "MISSENSE; G->V",
#///                              "NONSENSE; Y->STOP",
#///                              "NA; intron end", 
#///                              "MISSENSE; T->A",
#///                              "NONSENSE:MISSENSE R->Q",
#///                              "MISSENSE; S->C",
#///                              "MISSENSE; S->N",
#///                             "MISSENSE; L->P)")
write.csv(AnnotatedVcf_new,"AnnotatedVcf_new.csv")
write.csv(AnnotatedVcf_COSM,"AnnotatedVcf_COSM.csv")
write.csv(AnnotatedVcf_novel,"AnnotatedVcf_novel.csv")
##################################################################
### C) Reading in FPKM data and transform it into matrix format
##################################################################

#if under mac
#-------------Get FPKM Data-------------------------------------------------------------
fpkm_files_path <-paste0("/Users/yah2014/Documents/Programs/R/FPKM/", #Create fpkm_files_path file contains folder names
                         SampleList$Sample_ID,"_CuffLinks/genes.fpkm_tracking") 
FPKM <-data.frame(x= str(0), y= integer(0)) #Generate empty FPKM dataframe
FPKM_temp  <- read.csv(fpkm_files_path[1], header=T, sep="\t") #Read data from 1st sample
FPKM_temp <- FPKM_temp[order(FPKM_temp$tracking_id),] #Reorder the gene name
FPKM<-FPKM_temp$FPKM #Fill up FPKM dataframe with tracking id

#-------Get all FPKM Data----------
for(i in 2:nrow(SampleList)){ #skip the first one, which is added already
    FPKM_temp  <- read.csv(fpkm_files_path[i], header=T, sep="\t")
    FPKM_temp <- FPKM_temp[order(FPKM_temp$tracking_id),] #Reorder the gene name
    FPKM<-cbind(FPKM,FPKM_temp$FPKM)            
}

#intend to double check the tracking_id identity inside for loop
#give up because of speed.
#    for(n in 1:25286){         #Toooooo slow!!!
#        if(FPKM_temp[n,"tracking_id"]== FPKM[n,"tracking_id"])
#            FPKM[n,"i"]<-FPKM_temp[n,"FPKM"]
#    }
#-------name column and row----------------
dim(FPKM)
colnames(FPKM) <-SampleList_condition$Sample_ID # or annotation$SimpleLabel
rownames(FPKM) <-FPKM_temp$tracking_id

#########################################################
### D) Customizing and plotting the heat map
#########################################################


setwd("/Users/yah2014/Dropbox/Public/Olivier/R/WTCHG/Mutation")
getwd()
list.files()
# set a margin
par(oma=c(15,3,3,3))  # all sides have 3 lines of space  
par(mar=c(2,2,2,2))  
# ---QC------------
# ----------boxplot------------------------------------------------------
boxplot(log10(FPKM + 1),xlab="all 75 samples",cex=0.1, cex.lab=0.1, ylab = "log (base 10) RPKM + 1",main="FPKM for all Samples",las=2)

#-----cluster -------

par(mar=c(2,2,2,2))  
lm <-log(FPKM+1)
dm<-as.dist(1-cor(lm))
hm<-hclust(dm,method = "average")
p<-plot(hm,col = "black",cex.lab=0.1,main="Cluster Dendrogram for FPKM")
#-----heatmap -----------------------------------------------------

my_palette1 <- colorRampPalette(c("antiquewhite2", "green", "blue"))(n = 100)
my_palette2 <- colorRampPalette(c("red", "green"))(n = 2)
my_palette3 <- colorRampPalette(c("yellow", "orange", "red","green","blue","purple"))(n = 6)

# creat gene annotation label 
#//////mut_gene <- AnnotatedVcf_COSM[,!(colnames(AnnotatedVcf_COSM) %in% c("CHROM","POS","gene","ID","REF","ALT","QUAL","FILTER","INFO","counts",
#//////                                                                   "ZH2_rnaWT2","ZH4_rnaWT4"))] 
mut_gene <- AnnotatedVcf_novel[,!(colnames(AnnotatedVcf_novel) %in% c("CHROM","POS","gene_POS","gene","ID","REF","ALT","QUAL","FILTER","INFO","counts"))] 
mut_gene_clean <- mut_gene[,!(colnames(mut_gene) %in% c("ZH2_rnaWT2","ZH4_rnaWT4"))]

#-----heatmap---------------------



par(oma=c(10,3,3,3))
par(mar=c(10,2,2,2)+1)  
Heatmap(mut_gene, 
        column_title = "Cluster of RNA expression vs DNA mutation genes in WM",
        show_row_names = TRUE,
        col=my_palette1, 
        row_title_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 10),
        clustering_distance_columns = "euclidean",
        cluster_rows = FALSE,
        cluster_columns =  hm,
        row_dend_width = unit(5, "cm"),
        column_dend_height = unit(30, "mm"),
        heatmap_legend_param = list(title = "mutation rate %")
) 
