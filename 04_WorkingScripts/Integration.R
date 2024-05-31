#!/usr/local/bin/Rscript

# title: "Integration.R"

# function: 
#   Takes the Seurat Object post-PCA and runs 
#   Integration on the Seurat Object if the 
#   IntegrationMethod parameter is not set to 
#   NULL. Otherwise, this script will be skipped

#############################################################################################################################################################################

##################
# LOAD LIBRARIES #
##################

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

##################################
# READ IN PARAMS AND DIRECTORIES #
##################################

args <- commandArgs(trailingOnly = TRUE)

# RDS file from QC
params.SeuratObject <- args[1]


# Integration Method: Options( CCA, RPCA, Harmony, FastMNN, NULL) where NULL is to not run
params.IntegrationMethod <- args[2]

# Project Name
params.ProjectName <- args[3]


################
# Read in .rds # 
################
MergedSO <- readRDS(params.SeuratObject)


###################
# Run Integration #
###################
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


