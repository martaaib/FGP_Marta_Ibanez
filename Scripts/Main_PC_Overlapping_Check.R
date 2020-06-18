## ---Packages--- ##
library(dplyr)
library(readxl)

## ---Loading data of pregnancy complications--- ##
# setwd("~/Desktop/Universitat/3_ANY/TFG/Pregnancy Complications") , setting working directory
PC <- read_excel("Main_Pregnancy_Complication_Overlap_Secretome.xlsx")

## ---Pregnancy complications genes--- ##
# Selecting data
SGA <- PC[1:20,1] # Small for gestational age genes
colnames(SGA) <- "Gene"
LGA <- PC[1:6,2] # Large for gestational age genes
colnames(LGA) <- "Gene"
PE <- PC[1:91,3] # Preeclampsia genes
colnames(PE) <- "Gene"
IUGR <- PC[1:27,4] # Intrauterine growth restriction genes
colnames(IUGR) <- "Gene"
GDM <- PC[1:46,5] # Gestational diabetes mellitus genes
colnames(GDM) <- "Gene"
PC_check <- unique(rbind(SGA, LGA, PE, IUGR, GDM))

## ---Overlaying data--- ##

## SGA ##
SGA_LGA <- merge(SGA, LGA, by = "Gene") # 1
SGA_PE <- merge(PE, SGA, by = "Gene") # 17 
SGA_IUGR <- merge(SGA, IUGR, by = "Gene") # 6
SGA_GDM <- merge(SGA, GDM, by = "Gene") # 10

## LGA ##
LGA_PE <- merge(LGA, PE, by = "Gene") # 3
LGA_IUGR <- merge(LGA, IUGR, by = "Gene") # 3
LGA_GDM <- merge(LGA, GDM, by = "Gene") # 1

## PE ##
PE_GDM <- merge(PE, GDM, by = "Gene") # 28
PE_IUGR <- merge(IUGR, PE, by = "Gene") # 22

# GDM and IUGR ##
GDM_IUGR <- merge(GDM, IUGR, by = "Gene") # 10



