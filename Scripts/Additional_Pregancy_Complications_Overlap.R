## Intersection additional pregnancy complications and main pregnancy complications ##
require(BiocManager)
library(dplyr)
library(readxl)
setwd("~/Desktop/Universitat/3_ANY/TFG/Pregnancy Complications")

## Load data
## ---- Main Pregnancy Complications ---- ##
PE <- read_excel("PE.xlsx", col_names = FALSE)
GDM <- read_excel("GDM.xlsx", col_names = FALSE)
SGA <- read_excel("SGA.xlsx", col_names = FALSE)
LGA <- read_excel("LGA.xlsx", col_names = FALSE)
IUGR <- read_excel("IUGR.xlsx", col_names = FALSE)
num <- unique(rbind(PE,GDM,SGA,LGA,IUGR))
num <- unique(num)
## ---- Additional Pregnancy Complications ---- ##
PL <- read_excel("Preterm Labor.xlsx", col_names = FALSE)
MIS <- read_excel("Miscarriage.xlsx", col_names = FALSE)
PC <- read_excel("Pregnancy Complications Data Genes.xlsx")
total_unique <- PC[,9]
colnames(total_unique) <- "...1"
## --- Intersection between pregnancy complications --- ##
PE_PL <- merge(PE, PL, by = "...1")
GDM_PL <- merge(GDM, PL, by = "...1")
SGA_PL <- merge(SGA, PL, by = "...1")
LGA_PL <- merge(LGA, PL, by = "...1")
IUGR_PL <- merge(IUGR, PL, by = "...1")
PE_MIS <- merge(PE, MIS, by = "...1")
GDM_MIS<- merge(GDM, MIS, by = "...1")
SGA_MIS<- merge(SGA, MIS, by = "...1")
LGA_MIS<- merge(LGA, MIS, by = "...1")
IUGR_MIS<- merge(IUGR, MIS, by = "...1")
MIS_PL <- merge(PL, MIS, by = "...1")

## --- Prepare results  --- ## 
PL_list <- unique(rbind(PE_PL, GDM_PL, SGA_PL, LGA_PL, MIS_PL)) ## full list of PL overlapped with the PC
MIS_list <- unique(rbind(PE_MIS, GDM_MIS, SGA_MIS, LGA_MIS, MIS_PL)) ## full list of MIS overlapped with the PC
PL_specific <- anti_join(PL, PL_list, by = "...1")
MIS_specific <- anti_join(MIS, MIS_list, by = "...1")
MIS_PL_PC <- merge(MIS_list, PL_list, by = "...1")
MIS_PL2 <- merge(MIS_PL_PC, num, by = "...1") ## shared by main pc and additional
specific_MP <- anti_join(MIS_PL_PC, MIS_PL2, by = "...1") ## Shared only by PL and MIS
mis_pc <- anti_join(MIS_list,MIS_PL, by = "...1") ## Shared PC and MIS
Overlapps <- unique(rbind(MIS_PL,mis_pc))
try <- anti_join(num, total_unique, by = "...1")


