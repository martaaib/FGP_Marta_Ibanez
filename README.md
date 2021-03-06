# Defining the placental secretome
By: Marta Ibáñez

Supervised by: Dr Amanda Sferruzzi-Perri, Dr Tina Napso and Dr Xiaohui Zhao.

ESCI-UPF and University of Cambridge (Department of Physiology, Development and Neuroscience). Sferruzzi-Perri lab.

## For more information, go to: https://github.com/CTR-BFX/2020-Napso_Sferruzi-Perri  (private until the paper is published)

# Repository content:
## Readme:
Information about the repository and folders content.

## Data:
It contains some of the data used for the project:
- Ensemble_MGI_NCBI_human_mouse_Homolog_april_2020.csv: the homologs list.

- Pregnancy_Complications_Datasets.xlsx: datasets used for pregnancy complications.

- GEO_Control_Human_D8_N34673_GeneName_List.csv: list of genes from the GEO control Human datasets.

- GEO_Control_Mouse_D3_N47747_GeneName_List.csv: list of genes from the GEO control Mouse datasets.

## Results:
- Genes_Overlapping_Secretome_Additional_Pregnancy_Complications.xlsx: list of genes found in the 319 potentially secreted protein list and in published data from the additional pregnancy complications (Preterm labour and Miscarriage)

- Main_Pregnancy_Complication_Overlap_Secretome.xlsx : list of genes found in the 319 potentially secreted protein list and in already published literature of the selected main pregnancy complications (PE, GDM, Small for gestational age, Large for gestational age and IUGR).

- Overlap_Trophoblast_Organoid_final_Secretome_list.xlsx : list of genes found in the human trophoblast organoid and our 319 potentially secreted protein list.

## Scripts
- Additional_Pregancy_Complications_Overlap.R : Overlaying the miscarriage and preterm labour list of genes between each other and the main pregnancy complications and the 319 potentially secreted protein list.

- GO_data_Vis_secreted_list319.R : Data visualisation with "ggplot2" R package for the GO enrichment analysis results (Biological process, Molecular function and Cellular component).

- Main_PC_Overlapping_Check.R : Checking the overlap between the main pregnancy complications (PE, GDM, LGA, SGA and IUGR).

- Overlap_Trophoblast_Organoid_final_Secretome_list_AND_GO_Enrichment_Analysis.R: Overlaying the human trophoblast organoid data with our 319 potentially secreted protein list and data visualisation with "ggplot2" for the GO enrichment analysis.

- Protenomic_DataSort_PubCheck.R: Script used for the main peptidome analysis. It includes the three methodologies data (Conditioned medium, Sorted cells and cultured trophoblast). 

- Workflow_Numbers_Check.R: Script used for the checking of the main workflow final numbers and overlaying the different methodologies at different stages. 
