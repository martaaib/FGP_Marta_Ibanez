#!/usr/local/bin/Rscript
# R 3.6.2
#---------------------------------------------------------------------------------
# mass spectrum protenomic data analysis
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/2020-Napso_Sferruzi-Perri
#
#
# Analysis Performed by Xiaohui Zhao
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#---------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------


message("+---------------------------------------------------------------------------------------")
message("+                   Install useful packages and setup working directory                 ")
message("+---------------------------------------------------------------------------------------")

suppressPackageStartupMessages({
  library("dplyr")
  library("methods")
  library("utils")
  library("ggplot2")
  library("cowplot")
  library("Seurat")
  library("Matrix")
  library("useful")
  library("reshape2")
  library("biomaRt")
  library("scran")
  library("scater")
  library("SingleCellExperiment")
  library("bigmemory")
  library("mltools")
  library("rhdf5")
  library("monocle3")
  library("recommenderlab")
  library("readxl")
  library("GEOquery")
  library("UniProt.ws")
  library("rgdal")
  library("AnnotationDbi")
  library("org.Mm.eg.db")
  library("readxl")
})
options(future.globals.maxSize = 4000 * 1024^2)

Base.dir <- "/Users/xz289/Documents/CTR_ans48_0003"
Data.dir <- paste0(Base.dir, "/Original_Data")
Project <- "CTR_ans48_0003"
Out.dir <- paste0(Base.dir, "/Tina_Paper_Data")

message("+-----            General settings for three data sets           ----------------+")

ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description'), mart = ensembl)          
head(ensEMBL2id)
subEnsembl <- unique(ensEMBL2id[,2:3])
colnames(subEnsembl) <- c("MouseGeneName", "description")

## Check each data sets by using Uniprot retrive to get the proper UniprotKB entry name 

Revised.names <- paste0(c("F8VQ40", "F8WJE0", "E9PW66", "ES1", "CTGF", "SYHC", "CYR61",
                          "NIBL1","GRP78","GPR98","GPX41","YB039", "CS043", "CQ062", "GPX42",
                          "SK2L2", "USMG5", "CN166", "NIBAN", "K1468", "WISP1",
                          "F6VW30", "E9PW66","A0A0R4J026"), "_MOUSE")
NewRevi.names <- paste0(c("LAMA1", "SAMH1", "NP1L1", "GAL3A", "CCN2", "HARS1", "CCN1", 
                          "NIBA2", "BIP", "AGRV2","GPX4","MTLN","TRIR", "CYBC1","GPX4",
                          "MTREX", "ATPMD", "RTRAF", "NIBA1", "RELCH", "CCN4",
                          "1433T", "NP1L1", "FST"), "_MOUSE")
Rev.ProteinID <- c("P19137","Q60710","E9PW66","Q9D172","P29268","Q61035","P18406","Q8R1F1","P20029","Q8VHN7",
                   "O70325","Q8BT35","Q9D735","Q3TYS2","O70325","Q9CZU3","Q78IK2","Q9CQE8","Q3UW53","Q148V7",
                   "O54775","P68254", "P28656","P47931")

message("+--------------------------------------------------------------------------------------+")
message("+ Cultured_Primary: Read in different Raw Proteomic Data sets                          +")
message("+--------------------------------------------------------------------------------------+") 

PData3 <- read.csv(paste0(Data.dir, "/Cultured_Primary_Trophoblast_Data.csv"), header=T) ## 2131 proteins
PData3.sel <- PData3[, c(3,4,17,18:23,25)]
PData3.sel$ProteinID <-  unlist(lapply(as.character(PData3.sel$Accession), 
                                       function(x)  strsplit(x, split="[|]")[[1]][1]))
PData3.sel$Accession.Number <- unlist(lapply(as.character(PData3.sel$Accession), 
                                             function(x)  strsplit(x, split="[|]")[[1]][2]))
PData3.sel$Accession.New <- PData3.sel$Accession.Number
PData3.sel$MouseGeneName <- unlist(lapply(unlist(lapply(as.character(PData3.sel$Description),
                                                        function(x) strsplit(x, split="GN=")[[1]][2])),function(x) strsplit(x, split=" ")[[1]][1]))

for(i in 1:length(Revised.names)){
  rev.ind <- which(PData3.sel$Accession.Number==Revised.names[i])
  print(rev.ind)
  if(length(rev.ind)!=0){
    PData3.sel$Accession.New[rev.ind] <- NewRevi.names[i]
    PData3.sel$ProteinID[rev.ind] <- Rev.ProteinID[i]
  }
}

PData3.sel.QC <- subset(PData3.sel, X.Unique >=2 & X.10lgP >=12.7) 
## 1534, Protein-1533(GPX41,GPX42), Genes 1531 (Gnas, Naca, Gpx4), (Peptides <2, 597)
## Q9D735 (Trir) and Q3TYS2 (Cybc1) which has no gene names in the data,
## but found the gene names in uniport.

miss.ind <- which(is.na(PData3.sel.QC$MouseGeneName))
miss.ind.protein <- PData3.sel.QC[miss.ind,]

## Quality control, peptide >=2, 10lgP
PData3.sel.QC$MouseGeneName[miss.ind] <- c("Trir", "Cybc1")

## save the PData3.sel.QC for overlap use later
Protein.dat3c <- PData3.sel.QC[,c(13,14,11,10,4:8,2,3,1,12)]
colnames(Protein.dat3c) <- c("Accession.New", "MouseGeneName", "ProteinID", "Description", "Sample1", "Sample2", 
                             "Sample3","Sample4", "Sample5", "10lgP", "UniquePeptide", "Accession.Number.ori","Accession.Number")

## filter the proteins which expressed less than 3 samples
QC.submat <- PData3.sel.QC[,4:8]
QC.submatind <- ifelse(QC.submat==0,0,1)
QC.rind <- rowSums(QC.submatind)
PData3.sel.QC.GN.S4 <- PData3.sel.QC[QC.rind>=4,] ## 1208 (1206 unique genes)

PData3.final <- PData3.sel.QC.GN.S4[,c(13,14,11,10,4:8,2,3,1,12)]

colnames(PData3.final) <- colnames(Protein.dat3c) 
write.csv(PData3.final, file = paste0(Data.dir, "/Cultured_Trophoblast_N1208_G1206_Gnas_Naca_Filtered_Step2_Data.csv"), row.names=F, quote=T)


test3<- PData3.final[PData3.final$MouseGeneName%in%ensEMBL2id$external_gene_name==F,]
length(unique(test.sub$Accession.New)) ## 67
write.csv(test3, file = paste0(Data.dir, "/Cultured_Trophoblast_N67_G67_NotEnsem_Filtered_Step2_Data.csv"), row.names=F, quote=T)

message("+--------------------------------------------------------------------------------------+")
message("+ Conditonal_Medium:Read in different Raw Proteomic Data sets                          +")
message("+--------------------------------------------------------------------------------------+")

PData1 <- read.csv(paste0(Data.dir, "/Conditional_Medium_Data.csv"), header = T)
##  1467 x 14
PData1$Accession.Number.ori <- PData1$Accession.Number
Accession.Number <- PData1$Accession.Number
Protein.dat1 <- PData1[grepl("*_MOUSE*", Accession.Number)==T, ]
## 1445
Protein.dat1$Accession.Number <- unlist(lapply(as.character(Protein.dat1$Accession.Number.ori), function(x) strsplit(x, split=" ")[[1]][1]))
spindex <- which(grepl("^sp*", Protein.dat1$Accession.Number)==T)
Protein.dat1$Accession.Number[spindex] <- unlist(lapply(as.character(Protein.dat1$Accession.Number[spindex]), 
                                                        function(x) strsplit(x, split="[|]")[[1]][3]))

## Another two protein accession number with -DECOY
decoyindex <- which(grepl("*DECOY", Protein.dat1$Accession.Number)==T)
Protein.dat1$Accession.Number[decoyindex] <- unlist(lapply(as.character(Protein.dat1$Accession.Number[decoyindex]), 
                                                           function(x) strsplit(x, split="[-]")[[1]][1]))

Protein.dat1$MouseGeneName <- unlist(lapply(unlist(lapply(as.character(Protein.dat1$Identified.Proteins..1467.),
                                                          function(x) strsplit(x, split="GN=")[[1]][2])),function(x) strsplit(x, split=" ")[[1]][1]))
## 1445
## F6VW30, E9PW66,A0A0R4J026,F6QL70, D3Z536,A0A140T8L3, deleted by UniProt, april, 2020
## Ywhaq, Nap1l1,Fst, Gm17669,Gm8225,Rpl7a-ps5
## 1433T_MOUSE, NP1L1_MOUSE, FST_MOUSE

Protein.dat1$Accession.New <- Protein.dat1$Accession.Number

for(i in 1:length(Revised.names)){
  rev.ind <- which(Protein.dat1$Accession.Number==Revised.names[i])
  print(rev.ind)
  if(length(rev.ind)!=0){
    Protein.dat1$Accession.New[rev.ind] <- NewRevi.names[i]
  }
}


Protein.dat1c <- Protein.dat1[,c(18,17,4,10:14,16,5)]
colnames(Protein.dat1c) <- c("Accession.New", "MouseGeneName",  "Description", "Sample1", "Sample2", 
                             "Sample3","Sample4", "Sample5", "Accession.Number.ori","Accession.Number")

## 1445,--1443(KPYM, FBLN1),1441 (Psg16, Pkm, Fbln1,Tpm1)

missGindex <- which(is.na(Protein.dat1c$MouseGeneName))
## The protein TEST2  has no gene name corresponding.

Protein.dat1.newc <- Protein.dat1[-missGindex,c(18,17,4,10:14,16,5)] ## Remove TEST2 and othe useless columns.
## 1444
colnames(Protein.dat1.newc) <- colnames(Protein.dat1c)


## Filter the 4 out of 5 samples
QC.submat1 <- Protein.dat1.newc[,4:8]
QC.submatind1 <- ifelse(QC.submat1==0,0,1)
QC.rind1 <- rowSums(QC.submatind1)
Protein.dat1.new <- Protein.dat1.newc[QC.rind1>=4,]
## KPYM_MOUSE, sp|P52480|KPYM_MOUSE and sp|P52480-2|KPYM_MOUSE; D3YXQ6_MOUSE and D0VY58_MOUSE with gene Psg16
## 924
## Uniprot Retrive to get the UniProt ID.

## check
check.genes <- c("Ywhaq", "Nap1l1", "Fst", "Gm17669","Gm8225","Rpl7a-ps5")
Protein.dat1.new[unlist(lapply(check.genes, function(x) which(Protein.dat1.new$MouseGeneName==x))),]

Protein.dat1.uni <- read.csv(paste0(Data.dir, "/Conditional_Medium_UniProt_N924_filtered.csv"), header=T)[,1:2]
Protein.dat1.final <- merge(Protein.dat1.new,Protein.dat1.uni, by="Accession.New",all.x=T)
Protein.dat1.final[unlist(lapply(check.genes, function(x) which(Protein.dat1.final$MouseGeneName==x))),]
Protein.dat1.final[,11] <- as.character(Protein.dat1.final[,11])
Protein.dat1.final[6,11] <- "P68254"
## 924-protein 923 KPYM Gene 922 Psg16,Pkm

write.csv(Protein.dat1.final, file = paste0(Data.dir, "/Conditional_Medium_D924_G922_Psg16_Pkm_Filtered_Step2_Data.csv"),row.names =F)

test1 <- Protein.dat1.final[Protein.dat1.final$MouseGeneName%in%ensEMBL2id$external_gene_name==F,]
write.csv(test1, file = paste0(Data.dir, "/Conditional_Medium_Trophoblast_N29_G29_NotEnsem_Filtered_Step2_Data.csv"), row.names=F, quote=T)

message("+--------------------------------------------------------------------------------------+")
message("+ Sorted Trophoblast cell:Read in different Raw Proteomic Data sets                    +")
message("+--------------------------------------------------------------------------------------+")

Protein.dat2 <- read_excel(paste0(Data.dir, "/Sorted_Trophblast_Cells_Data.xlsx"))
colnames(Protein.dat2) <- c("Accession.Number", "S3", "S4", "S5", "P3", "P4", "P5")


Revised.names <- paste0(c("F8VQ40", "F8WJE0", "E9PW66", "ES1", "CTGF", "SYHC", "CYR61", "NIBL1","GRP78","GRP98","GPX41","YB039"), "_MOUSE")
NewRevi.names <- paste0(c("LAMA1", "SAMH1", "NP1L1", "GAL3A", "CCN2", "HARS1", "CCN1", "NIBA2", "BIP", "AGRV2","GPX4","MTLN"), "_MOUSE")

for(i in 1:length(Revised.names)){
  rev.ind <- which(Protein.dat2$Accession.Number==Revised.names[i])
  print(rev.ind)
  Protein.dat2$Accession.Number[rev.ind] <- NewRevi.names[i]
  Protein.dat2
}

Protein.dat2c <- Protein.dat2
## 1142

## Filter the proteins at least exist in 2 samples
QC.mat <- Protein.dat2c[,2:4]+Protein.dat2c[,5:7]
QC.mat.ind <- ifelse(QC.mat==0,0,1)
QC.mat.ind1 <- which(rowSums(QC.mat.ind)>=2)

Protein.dat2c.new <- Protein.dat2c[QC.mat.ind1,]
Protein.dat2c.new$Accession.New <- Protein.dat2c.new$Accession.Number
## Go to the Uniprot revive website to convert the entryname to Uniprot ID and GeneName
## Connect with the MouseGeneName using UniProt.ws library.
Protein.dat2.uni <- read.csv(paste0(Data.dir, "/Sorted_Trophoblast_UniProt_N682_filtered.csv"), header=T)[,1:4]

Protein.dat2.final <- merge(Protein.dat2c.new, Protein.dat2.uni, by = "Accession.New")
Protein.dat2.final$MouseGeneName <- as.character(Protein.dat2.final$MouseGeneName)
Protein.dat2.final[which(Protein.dat2.final$Accession.New=="H2B1K_MOUSE"),11] <- "Hist1h2bk"
colnames(Protein.dat2.final) <- c("Accession.New", "Accession.Number.ori", "S3", "S4", "S5", "P3", "P4", "P5", 
                                  "Accession.Number.uni", "ProteinID", "MouseGeneName")

write.csv(Protein.dat2.final, file = paste0(Data.dir, "/Sorted_Trophoblast_N682_G681_Naca_Filtered_Step2_Data.csv"),row.names =F)

test2 <- Protein.dat2.final[Protein.dat2.final$MouseGeneName%in%ensEMBL2id$external_gene_name==F,]
write.csv(test2, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N41_G41_NotEnsem_Filtered_Step2_Data.csv"), row.names=F, quote=T)


message("+----                 Perform Overlap with the GEO public Data           --------+")

MousePublic <- read.csv(paste0(Data.dir, "/GEO_Control_Mouse_D3_N47747_GeneName_List.csv"), header=T)
colnames(MousePublic) <- "MouseGeneName"

Dat1.overlap1 <- Protein.dat1.final[Protein.dat1.final$MouseGeneName%in%MousePublic[,1]==T,] 
## 908-907--906, Psg16, Pkm
Dat2.overlap1 <- Protein.dat2.final[Protein.dat2.final$MouseGeneName%in%MousePublic[,1]==T,] 
## 642-642--641, Naca
Dat3.overlap1 <- PData3.final[PData3.final$MouseGeneName%in%MousePublic[,1]==T,]
## 1180-1180-1178, Gnas, Naca, 28 genes are not found in ensEMBLE and also not found in public Mouse.


message("+----                 Perform Homolog overlap           --------+")

HumanPublic <- read.csv(paste0(Data.dir, "/GEO_Control_Human_D8_N34673_GeneName_List.csv"), header=T)
HomoData <- read.csv(paste0(Data.dir, "/", Project, "-Ensemble_MGI_NCBI_human_mouse_Homolog_april_2020.csv"), header=T)

HomoData.overPub <- HomoData[as.character(HomoData$HumanGeneName)%in%as.character(HumanPublic[,1])==T,]
## make HomoData.overPub unique for mouse and Human matching and remove all of the 1-more or more-1 genes
Mmulit <- HomoData.overPub[duplicated(HomoData.overPub$MouseGeneName)==T,]
## 5362
Hmulit <- HomoData.overPub[duplicated(HomoData.overPub$HumanGeneName)==T,]
## 6679

HomoData.overPub.Muni <- unique(as.character(HomoData.overPub$MouseGeneName)) 
colnames(HomoData.overPub.Muni) <- "MouseGeneName"

Dat1.overlap2 <- Dat1.overlap1[Dat1.overlap1$MouseGeneName%in%HomoData.overPub.Muni==T,] 
## 898-898--896, Psg16, Pkm
Dat2.overlap2 <- Dat2.overlap1[Dat2.overlap1$MouseGeneName%in%HomoData.overPub.Muni==T,] 
## 638-638--637, Naca
Dat3.overlap2 <- Dat3.overlap1[Dat3.overlap1$MouseGeneName%in%HomoData.overPub.Muni==T,] 
## 1175-1175-1173, Gnas, Naca

message("+----  Suggestions from Tina to remove some genes, due to multiple homologus with human       ----+")

Filter.proteins <- read.csv(paste0(Data.dir, "/Filter_Proteins_List.csv"), header=T)

Dat1.overlap2.new <- Dat1.overlap2[Dat1.overlap2$MouseGeneName%in%Filter.proteins$MouseGeneName==F,] 
## 881-881--879, Psg16, Pkm
Dat2.overlap2.new <- Dat2.overlap2[Dat2.overlap2$MouseGeneName%in%Filter.proteins$MouseGeneName==F,] 
## 634-634--633, Naca
Dat3.overlap2.new <- Dat3.overlap2[Dat3.overlap2$MouseGeneName%in%Filter.proteins$MouseGeneName==F,] 
## 1172-1172-1170, Gnas, Naca

test1.mh.over <- Dat1.overlap1[Dat1.overlap1$MouseGeneName%in%Dat1.overlap2.new$MouseGeneName==T,]
test1.mh.nover <- Dat1.overlap1[Dat1.overlap1$MouseGeneName%in%Dat1.overlap2.new$MouseGeneName==F,]

test2.mh.over <- Dat2.overlap1[Dat2.overlap1$MouseGeneName%in%Dat2.overlap2.new$MouseGeneName==T,]
test2.mh.nover <- Dat2.overlap1[Dat2.overlap1$MouseGeneName%in%Dat2.overlap2.new$MouseGeneName==F,]

test3.mh.over <- Dat3.overlap1[Dat3.overlap1$MouseGeneName%in%Dat3.overlap2.new$MouseGeneName==T,]
test3.mh.nover <- Dat3.overlap1[Dat3.overlap1$MouseGeneName%in%Dat3.overlap2.new$MouseGeneName==F,]

message("+---write out each step proteins list for flowchart Data1, 2,3 ----------------------+")

write.csv(PData3.sel.QC, file = paste0(Data.dir, "/Cultured_Trophoblast_N1534_G1531_Gnas_Naca_Gpx4_Filtered_pep2_Step1_Data.csv"), row.names=F, quote=T)
write.csv(Dat3.overlap1, file = paste0(Data.dir, "/Cultured_Trophoblast_N1180_G1178_Gnas_Naca_Filtered_pep2_Step3_MousePub_Data.csv"), row.names=F, quote=T)
write.csv(Dat3.overlap2.new, file = paste0(Data.dir, "/Cultured_Trophoblast_N1172_G1170_Gnas_Naca_Filtered_pep2_Step3_HumPub_Data.csv"), row.names=F, quote=T)
write.csv(test3.mh.nover,file = paste0(Data.dir, "/Cultured_Trophoblast_N8_G8_Filtered_pep2_Step4_MusUni_Data.csv"), row.names=F, quote=T )

write.csv(Protein.dat2c, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N1142_Step1_Data.csv"), row.names=F, quote=T)
write.csv(Dat2.overlap1, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N642_G641_Naca_Filtered_Step3_MousePub_Data.csv"), row.names=F, quote=T)
write.csv(Dat2.overlap2.new, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N634_G633_Naca_Filtered_Step3_HumPub_Data.csv"), row.names=F, quote=T)
write.csv(test2.mh.nover, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N8_G8_Filtered_Step4_MusUni_Data.csv"), row.names=F, quote=T)


write.csv(Protein.dat1c, file = paste0(Data.dir, "/Conditioned_Medium_Trophoblast_N1445_G1441_Psg16_Pkm_Tpm1_Fbln1_Step1_Data.csv"), row.names=F, quote=T)
write.csv(Dat1.overlap1, file = paste0(Data.dir, "/Conditioned_Medium_Trophoblast_N881_G879_Psg16_Pkm_Filtered_Step3_MousePub_Data.csv"), row.names=F, quote=T)
write.csv(test1.mh.nover, file = paste0(Data.dir, "/Conditioned_Medium_Trophoblast_N27_G27_Filtered_Step4_MusUni_Data.csv"), row.names=F, quote=T)

message("+-------------                   END                 --------------------------+")