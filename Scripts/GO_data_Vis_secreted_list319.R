## GO enrichment Data Visualisation 320 secreted genes list
# setwd("~/Desktop/Universitat/3_ANY/Lab/Lab_meeting")
library(ggplot2)

## Biological Process
Bio_Process <- read.delim("Bio_Process_319_list.txt")
colnames(Bio_Process)[3] <- "Number.of.Genes"
colnames(Bio_Process)[7] <- "P.value"
Bio_Process <- Bio_Process[Bio_Process[7] < 0.05,]
Bio_Process <- Bio_Process[order(Bio_Process$Number.of.Genes, decreasing = TRUE),]
ggplot(data = Bio_Process[4:14,], mapping = aes(x=reorder(GO.biological.process.complete, Number.of.Genes), y = as.factor(Number.of.Genes))) + geom_col( fill = "darkseagreen3") + coord_flip() + theme_minimal(base_size = 15) +
  xlab("Biological Process") + ylab("Number of genes") + geom_text(aes(label = Number.of.Genes), position = "stack") + ggtitle("GO terms Secreted Data") 

## Molecular Function
Mol_Function <- read.delim("Mol_Func_319_List.txt")
colnames(Mol_Function)[3] <- "Number.of.Genes"
colnames(Mol_Function)[7] <- "P.value"
Mol_Function <- Mol_Function[Mol_Function[7] < 0.05,]
Mol_Function <- Mol_Function[order(Mol_Function$Number.of.Genes, decreasing = TRUE),]
ggplot(data = Mol_Function[2:12,], mapping = aes(x=reorder(GO.molecular.function.complete, Number.of.Genes), y = as.factor(Number.of.Genes))) + geom_col( fill = "aquamarine3") + coord_flip() + theme_minimal(base_size = 15) +
  xlab("Molecular Function") + ylab("Number of genes") + geom_text(aes(label = Number.of.Genes), position = "stack") + ggtitle("GO terms Secreted Data")

## Cellular Component
Cel_Comp <- read.delim("Cell_Comp_319_list.txt")
colnames(Cel_Comp)[3] <- "Number.of.Genes"
colnames(Cel_Comp)[7] <- "P.value"
Cel_Comp <- Cel_Comp[Cel_Comp[7] < 0.05,]
Cel_Comp <- Cel_Comp[order(Cel_Comp$Number.of.Genes, decreasing = TRUE),]
Cel_Comp <- Cel_Comp[-c(2),]
ggplot(data = Cel_Comp[2:12,], mapping = aes(x=reorder(GO.cellular.component.complete, Number.of.Genes), y = as.factor(Number.of.Genes))) + geom_col( fill = "cadetblue3") + coord_flip() + theme_minimal(base_size = 15) +
  xlab("Cellular Component") + ylab("Number of genes") + geom_text(aes(label = Number.of.Genes), position = "stack") + ggtitle("GO terms Combined Secreted Data")
