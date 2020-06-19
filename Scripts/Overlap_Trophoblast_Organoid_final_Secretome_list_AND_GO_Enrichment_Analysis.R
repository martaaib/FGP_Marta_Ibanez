## Trophoblast Organoid Secretome Data ##
setwd("~/Desktop/Universitat/3_ANY/Lab/Trophoblast_new_secreted_data")
library(VennDiagram)
library(readxl)
library(xlsx)
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
library(grid)
## Load T.O data
trophoblast_org_data <- read_excel("homologs_edited.xlsx")
troph_org <- trophoblast_org_data[,c(1,3)]
colnames(troph_org) <- c("Gene Name", "Protein Name")

## New secretome
secretome <- read_excel("320-final list.xlsx")
secretome <- apply(secretome,2,toupper) ## converting into upper case letters
secretome <- data.frame(secretome)
colnames(secretome) <- "Gene Name"

## Merge
TO_secretome <- merge(troph_org, secretome, by = "Gene Name")
## Write excel file
write.xlsx2(TO_secretome,"Overlap_Trophoblast_Organoid_Secretome.xlsx")



## Venn Diagram ##
V <- venn.diagram(
  x = c(secretome, troph_org[,1]),
  category.names = c("Secretome" , "Trophoblast Human Organoid"),
  filename = NULL,
  #output = TRUE ,
  #imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 1200,
  compression = "lzw",
  lwd = 1,
  col=c("paleturquoise3", 'royalblue3'),
  fill = c(alpha("paleturquoise3",0.3), alpha('royalblue3',0.3)),
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("paleturquoise3", 'royalblue3'),
  # Numbers
  cex = .8,
  fontface = "bold",
  fontfamily = "sans",
)
grid.newpage()
grid.draw(V)

## Data Visualisation on GO enrichment analysis ##
Bio_process <- read.delim("Bio_process_TO-Sec.txt")
colnames(Bio_process)[3] <- "Number.of.Genes"
colnames(Bio_process)[7] <- "P.value"
Bio_process<- Bio_process[order(-Bio_process[,3], Bio_process[,7]), ]
ggplot(data = Bio_process[1:10,], mapping = aes(x=reorder(GO.biological.process.complete, Number.of.Genes), y = as.factor(Number.of.Genes))) + geom_col( fill = "steelblue3") + coord_flip() + theme_minimal(base_size = 20) +
  xlab("Biological Process") + ylab("Number of genes") + geom_text(aes(label = Number.of.Genes), position = "stack") + ggtitle("GO terms HTO and Secreted Data") 

Mol_Function <- read.delim("Mol_Func_TO-Sec.txt")
colnames(Mol_Function)[3] <- "Number.of.Genes"
colnames(Mol_Function)[7] <- "P.value"
Mol_Function<- Mol_Function[order(-Mol_Function[,3], Mol_Function[,7]), ]
ggplot(data = Mol_Function[2:11,], mapping = aes(x=reorder(GO.molecular.function.complete, Number.of.Genes), y = as.factor(Number.of.Genes))) + geom_col( fill = "lightcyan3") + coord_flip() + theme_minimal(base_size = 20) +
  xlab("Molecular Function") + ylab("Number of genes") + geom_text(aes(label = Number.of.Genes), position = "stack") + ggtitle("GO terms HTO and Secreted Data") 


Cel_Comp <- read.delim("Overlapped_Cel_Comp.txt")
colnames(Cel_Comp)[3] <- "Number.of.Genes"
colnames(Cel_Comp)[7] <- "P.value"
Cel_Comp<- Cel_Comp[order(-Cel_Comp[,3], Cel_Comp[,7]), ]
ggplot(data = Cel_Comp[2:11,], mapping = aes(x=reorder(GO.cellular.component.complete, Number.of.Genes), y = as.factor(Number.of.Genes))) + geom_col( fill = "thistle3") + coord_flip() + theme_minimal(base_size = 20) +
  xlab("Cellular Component") + ylab("Number of genes") + geom_text(aes(label = Number.of.Genes), position = "stack") + ggtitle("GO terms HTO and Secreted Data") 








