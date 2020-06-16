## Final checkig number of diagrams
setwd("~/Desktop/Universitat/3_ANY/Lab/Data3/Finish_Workflow_diagram_numbers")

## Checking data1 ##
data1 <- read.csv("CTR_ans48_0003-Cultured_Mouse_Human_Public_common_N873_G872_Psg16_FilterPrl_24_03_20.csv")
nrow(unique(data.frame(data1[,1]))) ## 872 genes
data1_genes <- unique(data.frame(data1[,1:2]))
colnames(data1_genes) <- "MouseGeneName"

## Data1 layer 1 Gene Level, 918 ##
data1_918 <- read.csv("CTR_ans48_0003-Cultured_F4out5_peptide2_N919_G918_Psg16_21_05_19.csv")
data1_918 <- data1_918[,8:9]
data1_918_genes <- unique(data.frame(data1_918[,2]))
colnames(data1_918_genes) <- "MouseGeneName"

## Data1 layer 2 Gene Level, 900 ##
data1_900 <- read.csv("CTR_ans48_0003-Cultured_Mouse_Public_common_N901_G900_Psg16_21_05_19.csv")
data1_900_genes <- unique(data.frame(data1_900[,2]))
colnames(data1_900_genes) <- "MouseGeneName"

## Checking data2
data2 <- read.csv("CTR_ans48_0003-Sorted_Mouse_Human_Public_common_N655_G654_Naca_FilterPrl_24_03_20.csv")
nrow(unique(data.frame(data2[,1]))) ## 654 genes
data2_genes <- unique(data.frame(data2[,1]))
colnames(data2_genes) <- "MouseGeneName"

## Data2 layer 1 Gene Level, 680 ##
data2_680 <- read.csv("CTR_ans48_0003-Isolated_F2out3_peptide2_N681_G680_Naca_21_05_19.csv")
data2_680 <- data2_680[,8:9]
data2_680_genes <- unique(data.frame(data2_680[,2]))
colnames(data2_680_genes) <- "MouseGeneName"

## Data2 layer 2 Gene Level, 662 ##
data2_662 <- read.csv("CTR_ans48_0003-Isolated_Mouse_Public_common_N663_G662_Naca_21_05_19.csv")
data2_662_genes <- unique(data.frame(data2_662[,2]))
colnames(data2_662_genes) <- "MouseGeneName"

## Checking data3 ##
data3 <- read.csv("CTR_ans48_0003-Data3_Mouse_Human_Public_common_N1169_G1167_Naca_Gnas_FilterPrl_24_03_20.csv")
nrow(unique(data.frame(data3[,1]))) ## 1167 genes
data3_genes <- unique(data.frame(data3[,1:2]))
colnames(data3_genes) <- "MouseGeneName"

## Data3 layer 1 Gene Level, 1206 ##
load("PData3.RData")
data3_1208 <- PData3.sel.QC.GN.S4[,11:12]
data3_1208_genes <- unique(data.frame(data3_1208[,2]))
data3_1206_genes <- data3_1208_genes
colnames(data3_1206_genes) <- "MouseGeneName"


## Data3 layer 2 Gene Level, 1177 ##

specific_data3 <- read.csv("Data3_Mouse_Human_public_Munique_NP6_NG6.csv")
specific_data3_genes <- specific_data3[,c(1,12)]


## Comparison data1-data2
data1_data2_genes <- merge(data1_genes, data2_genes, by = "MouseGeneName")
write.csv(data1_data2_genes, "data1N872_data2N654_CommonN362.csv")
## Comparison data1-data3
data1_data3_genes <- merge(data1_genes, data3_genes, by = c("MouseGeneName"))
write.csv(data1_data3_genes, "data1N872_data3N166_CommonN545.csv")
## Comparison data2-data3
data2_data3_genes <- merge(data2_genes, data3_genes, by = "MouseGeneName")
write.csv(data2_data3_genes, "data2N654_data3N1667_CommonN519.csv")

## Last layer ##
data1_data2_data3_commongenes <- merge(data1_data2_genes, data3_genes, by = "MouseGeneName")
write.csv(data1_data2_data3_commongenes, "data1N872_data2N654_data3N1667_CommonN334.csv")

## Comparison data1_918- data2_680
data1_918_data2_680 <- merge(data1_918_genes, data2_680_genes, by = "MouseGeneName")
write.csv(data1_918_data2_680, "data1N918_data2N680_CommonN379.csv")

## Comparison data1_918- data3_1206
data1_918_data3_1206 <- merge(data1_918_genes, data3_1206_genes, by = "MouseGeneName")
write.csv(data1_918_data3_1206, "data1N918_data3N1206_CommonN563.csv")

## Comparison data2_680- data3_1206
data2_680_data3_1206 <- merge(data2_680_genes, data3_1206_genes, by = "MouseGeneName")
write.csv(data2_680_data3_1206, "data2N680_data3N1206_CommonN537.csv")

## First layer ##
data1_2_3_first_layer_common <- merge(data1_918_data2_680, data3_1206_genes, by = "MouseGeneName")
write.csv(data1_2_3_first_layer_common, "data1N918_data2N680_data3N1206_CommonN345.csv")

## Comparison data1_900- data2_662
data2_662_data1_900 <- merge(data2_662_genes, data1_900_genes, by = "MouseGeneName")
write.csv(data2_662_data1_900, "data1N900_data2N662_CommonN367.csv")

## Comparison data1_900- data3_1177 + data3_1169
data3_1169_data1_900 <- merge(data1_900_genes, data3_genes, by = "MouseGeneName")
data1_900_specific3 <- merge(data1_900_genes, specific_data3, by = "MouseGeneName")
prova <- merge(data1_900_specific3, data3_1169_data1_900, by = "MouseGeneName", all.x = TRUE, all.y = TRUE)
write.csv(prova, "data1N900_data3N1173_CommonN549.csv")

## Comparison data2_662- data3_1177 + data3_1169
data3_1169_data2_662 <- merge(data2_662_genes, data3_genes, by = "MouseGeneName")
data2_662_specific3 <- merge(data2_662_genes, specific_data3, by = "MouseGeneName")
prova2 <- merge(data3_1169_data2_662, data2_662_specific3, by = "MouseGeneName", all.x = TRUE, all.y = TRUE)
write.csv(prova2, "data2N662_data3N1173_CommonN522.csv")

## Second Layer ##
data123_sec <- merge(data2_662_data1_900, data3_genes, by = "MouseGeneName" )
prova3 <- merge(data2_662_data1_900, specific_data3, by = "MouseGeneName")
write.csv(data123_sec, "data1N900_data2N662_data3_1173_CommonN336.csv")

## Protein Level ##

## 1 and 2 ##
data1_data2_prots <- merge(data1_918, data2_680, by = "MouseGeneName")
write.csv(data1_data2_prots, "data1N919_data2N681_CommonN380.csv")

## 1 and 3 ##
data1_data3_prots <- merge(data1_918, data3_1208, by = "MouseGeneName")
write.csv(data1_data3_prots, "data1N919_data3N1208_CommonN564.csv")

## 2 and 3 ##
data2_data3 <- merge(data2_680, data3_1208, by = "MouseGeneName")
write.csv(data2_data3, "data2N681_data3N1208_CommonN541.csv")

## Protein level all lists ##
data1_2_3_prots <- merge(data1_data2_prots, data3_1208, by = "MouseGeneName")
write.csv(data1_2_3_prots, "data1N919_data2N681_data3N1208_CommonN348.csv")

