#!/usr/local/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                        Title: Merge.R                                      ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes all of the Seurat Objects post DoubletFinder.R that belong to a specific         ═╣
# ╠═     analysis group and run them through the following steps: Merge, Find Variable Features ═╣
# ╠═     Generate merged Quality Control Plots, and Scaling.                                    ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝
# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(dplyr)
library(Matrix)
library(viridis)
library(tidyverse)
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(SeuratWrappers)
library(Seurat.utils)
library(SingleCellExperiment)
library(gprofiler2)
library(stringr)
library(patchwork)
library(presto)

# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
args <- commandArgs(trailingOnly = TRUE)

# Project Name for this analysis group
params.ProjectName <- args[1]

# Variable To Regress For Scaling
params.VarsToRegress <- args[2]
params.VarsToRegress <- unlist(strsplit(params.VarsToRegress, ","))
print(params.VarsToRegress)

# ╔════════════════════╗
# ╠═ Load Seurat .rds ═╣
# ╚════════════════════╝
SeuratRDSfiles <- list.files(pattern = "\\.rds$", ignore.case = T)

SampleNames <- c()

counter <- 1 
for( file in SeuratRDSfiles){
    ObjName <- gsub("_.*","", file)
    SampleNames <- append(SampleNames,ObjName)
    assign(ObjName, readRDS(file))
}

count <- 2
countMax <- length(SampleNames)
SOlist2up <- c()
while (count <= countMax){
    SOlist2up <- append(SOlist2up, get(SampleNames[count]))
    count <- count + 1
}

MergedSO <- merge(get(SampleNames[1]), y= SOlist2up, add.cell.ids = SampleNames, project = params.ProjectName)
MergedSO

# ╔═══════════════════════════════════════════╗
# ╠═ Normalize Data & Find Variable Features ═╣
# ╚═══════════════════════════════════════════╝
MergedSO <- NormalizeData(MergedSO)
MergedSO <- FindVariableFeatures(MergedSO)

# ╔═════════════════════════════════════╗
# ╠═ Generate Merged Post-Filter Plots ═╣
# ╚═════════════════════════════════════╝
plot1 <- FeatureScatter(MergedSO, feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle = T, group.by = "orig.ident")
plot2 <- FeatureScatter(MergedSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle = T, group.by = "orig.ident")
vplot1 <- VlnPlot(MergedSO, features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1, group.by = "orig.ident")

pdf(file = paste0(params.ProjectName,'postMergeQC.pdf'), title = paste0(params.ProjectName,' Post Merge QC'), width = 11, height = 8.5)
plot1 + plot2
vplot1
dev.off()

# ╔══════════════════════════════╗
# ╠═ Scale Merged Seurat Object ═╣
# ╚══════════════════════════════╝
all.genes <- rownames(MergedSO)
MergedSO = ScaleData(MergedSO, features = all.genes ,vars.to.regress = params.VarsToRegress)

# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
saveRDS(MergedSO, paste0(params.ProjectName,"Merged_SO.rds"))

# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
sink(paste0(params.ProjectName,"_Mergedvalidation.log"))
print(MergedSO)
sink()