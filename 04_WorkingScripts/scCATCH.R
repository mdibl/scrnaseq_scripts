#!/usr/local/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                      Title: scCATCH.R                                      ═╣
# ╠═                                   Updated: August 1 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes in the Seurat Object generated from FindNeighborsClustersMarkers.R. Uses         ═╣
# ╠═     scCATCH to identify cell identities for a specified cell clustering.                   ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝
# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(dplyr)
library(stringr)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(scCATCH)

# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
args <- commandArgs(trailingOnly = TRUE)

# RDS file from QC
params.SeuratObject <- args[1]

# Resolutions
params.Resolution <- args[2]
if (toupper(params.Resolution) == "NULL") {
    params.Resolution <- "0.3"
}

# Integration Method: Options( CCA, RPCA, Harmony, FastMNN, NULL) where NULL is to not run
params.IntegrationMethod <- args[3]
if (toupper(params.IntegrationMethod) == "NULL"){
    params.IntegrationMethod <- "unintegrated"
}

# Organism Needs to match DB
params.Organism <- args[4]

# Tissue List https://github.com/ZJUFanLab/scCATCH/wiki
params.Tissue <- args[5]
params.Tissue <- as.character(unlist(strsplit(params.Tissue, ",")))
print(params.Tissue)

params.scaleMethod <- args[6]

# Project Name
params.ProjectName <- args[7]

# Custom scCATCH file
#params.Database <- args[4]

# ╔══════════════════╗
# ╠═ Organism Check ═╣
# ╚══════════════════╝
message1 <- NULL
if (toupper(params.Organism) == "HUMAN" | 
    toupper(params.Organism) == "HSAPIENS" | 
    toupper(params.Organism) == "H SAPIENS" | 
    toupper(params.Organism) == "HOMO SAPIENS" | 
    toupper(params.Organism) == "HOMOSAPIENS" |
    toupper(params.Organism) == "HS") {
    message1 <- "NOTE:SCCATCH:Running with Human"
    message(message1)
    params.Organism <- "Human"
}else if (toupper(params.Organism) == "MOUSE" | 
          toupper(params.Organism) == "MMUSCULUS" | 
          toupper(params.Organism) == "M MUSCULUS" | 
          toupper(params.Organism) == "MUS MUSCULUS" | 
          toupper(params.Organism) == "MUSMUSCULUS" |
          toupper(params.Organism) == "MM") {
    message1 <- "NOTE:SCCATCH:Running with Mouse"
    message(message1)
    params.Organism <- "Mouse"
}else {
    message1 <- "ERROR:SCCATCH:Cell Identity Prediction only compatible with Human and Mouse"
    message(message1)
    quit(status = 1)    
}

# ╔════════════════╗
# ╠═ Tissue Check ═╣
# ╚════════════════╝
message2 <- NULL
cellMatchOrg <- subset(cellmatch, subset = cellmatch$species == params.Organism)

OrgTissues <- unique(sort(cellMatchOrg$tissue))
for (i in params.Tissue ){
    if (i %in% OrgTissues == F){
        message2 <- paste0("ERROR:SCCATCH:The following tissue is not in the database for ", params.Organism,": ", i)
        message(message2)
        quit(status = 1)    
    }
}

# ╔════════════════════╗
# ╠═ Load Seurat .rds ═╣
# ╚════════════════════╝
MergedSO <- readRDS(params.SeuratObject)

# ╔═════════════════════════╗
# ╠═ Create scCATCH Object ═╣
# ╚═════════════════════════╝
if (params.scaleMethod == "SD"){ 
    MergedDGC <- MergedSO[['RNA']]$data
}else {
    MergedDGC <- MergedSO[['SCT']]$data
}
MergedDGC@Dimnames$features <- rownames(MergedSO)
MergedDGC@Dimnames$cells <- colnames(MergedSO)
MergedDGC@Dimnames <- MergedDGC@Dimnames[-c(1,2)]

# ╔══════════════════════╗
# ╠═ Run scCATCH Object ═╣
# ╚══════════════════════╝
MergedDGC <- rev_gene(data = MergedDGC, data_type = "data", species = params.Organism, geneinfo = geneinfo)
Clusterings <- as.character(MergedSO[[paste0(params.IntegrationMethod,"Res.",params.Resolution)]][,1])
MergedSCCO <- createscCATCH(data = MergedDGC, cluster = Clusterings)
MergedSCCO <- findmarkergene(object = MergedSCCO, species = params.Organism, marker = cellmatch, tissue = params.Tissue)

MergedSCCO <- findcelltype(object = MergedSCCO)

# ╔══════════════════════════════════════╗
# ╠═ Add Cell Identity to Seurat Object ═╣
# ╚══════════════════════════════════════╝
MergedSO$CIscCATCH <- MergedSO[[paste0(params.IntegrationMethod,"Res.",params.Resolution)]][,1]

Cellframe <- as.data.frame(list(Cluster = MergedSO$CIscCATCH))
Cellframe$Cluster <- as.character(Cellframe$Cluster)
Cellframe$RowNum <- 1:nrow(Cellframe) 
Cellframe$RowNum <- str_pad(Cellframe$RowNum, nchar(nrow(Cellframe)), pad = "0" )
CIframe <- as.data.frame(list(Cluster = as.integer(MergedSCCO@celltype$cluster), CellType = MergedSCCO@celltype$cell_type), row.names = NULL)
Namedframe <- merge(Cellframe,CIframe, by.x="Cluster", by.y="Cluster", all.x = T)
Namedframe$CellType <- ifelse(is.na(Namedframe$CellType), Namedframe$Cluster, Namedframe$CellType)
Namedframe <- Namedframe[order(Namedframe$RowNum),]

MergedSO$CIscCATCH <- Namedframe$CellType

# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
SaveSeuratRds(MergedSO, file = paste0("07",params.ProjectName, "_scCATCHSO.rds"))

# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
sink(paste0("07_",params.ProjectName,"_scCATCHValidation.log"))
print("╔══════════════════════════════════════════════════════════════════════════════════════════════╗")
print("╠  scCATCH.R log")
print(paste0("╠  Analysis Group: ", params.ProjectName))
print("╚══════════════════════════════════════════════════════════════════════════════════════════════╝")
print("Seurat Object Status:")
print(MergedSO)
print("Cell Identity Predictions:")
print(CIframe)
cat("\n")
if(length(message1)!= 0){
    print(message1)
}
if(length(message2)!= 0){
    print(message2)
}
sink()

sink(paste0("07_",params.ProjectName,"_scCATCHVersions.log"))
print("╔══════════════════════════════════════════════════════════════════════════════════════════════╗")
print("╠  scCATCH.R Versions")
print(paste0("╠  Analysis Group: ", params.ProjectName))
print("╚══════════════════════════════════════════════════════════════════════════════════════════════╝")
sessionInfo()
sink()
