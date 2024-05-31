#!/usr/local/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                     Title: Integration.R                                   ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes in the Seurat Object generated from RunPCA.R. Runs integration on the Seurat     ═╣
# ╠═     Object if the Integration Method parameter is not set to NULL. Otherwise, this script  ═╣
# ╠═     is skipped.                                                                            ═╣
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
library(ggplot2)
library(ggsankey)
library(DoubletFinder)
library(stringr)
library(patchwork)
library(loupeR)
library(presto)

# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
args <- commandArgs(trailingOnly = TRUE)

# RDS file from QC
params.SeuratObject <- args[1]

# Integration Method: Options( CCA, RPCA, Harmony, FastMNN, NULL) where NULL is to not run
params.IntegrationMethod <- args[2]

# Project Name
params.ProjectName <- args[3]

# ╔════════════════════╗
# ╠═ Load Seurat .rds ═╣
# ╚════════════════════╝
MergedSO <- readRDS(params.SeuratObject)

# ╔═══════════════════╗
# ╠═ Run Integration ═╣
# ╚═══════════════════╝
if(params.IntegrationMethod == "FastMNN"){
    MergedSO <- IntegrateLayers(object = MergedSO, method = paste0(params.IntegrationMethod,"Integration"), new.reduction = paste0("integrated.",params.IntegrationMethod))
}else{
    MergedSO <- IntegrateLayers(object = MergedSO, method = paste0(params.IntegrationMethod,"Integration"), orig.reduction = "pca", new.reduction = paste0("integrated.",params.IntegrationMethod))
}

######################
# Save Seuart Object #
######################
SaveSeuratRds(MergedSO, file = paste0(params.ProjectName, "_Integrated.rds"))

sink(paste0(params.ProjectName,"validation.log"))

MergedSO

sink()


