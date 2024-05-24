#!/usr/local/bin/Rscript

# title: "Doublet_Finder.R"

# function: 
#   Takes the single sample Seurat objects through a 
#   process of finding doublets (Using Doublet Finder)
#   and removing the predicted doublets based on the 
#   expected number of cells for the sample. 


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

# Variables To Regress in the Scaling
params.VarsToRegress <- args[2]
params.VarsToRegress <- unlist(strsplit(params.VarsToRegress, ","))
print(params.VarsToRegress)

# Number of PCs to use, default to auto
params.PCs <- args[3]
if (params.PCs != "auto"){
    params.PCs <- as.integer(params.PCs)
}

# Data Dir
params.DataDir <- args[4]

# SampleName
params.sampleName <- args[5]

#############################
# Find Doublets Preperation #
#############################

# read in .rds  
samp <- readRDS(params.SeuratObject)

# create a list of all genes
all.genes <- rownames(samp)

# Scale Data
samp <- ScaleData(samp, features = all.genes ,vars.to.regress = params.VarsToRegress)


# Run PCA on individual samples
samp <- RunPCA(samp, features = all.genes ,npcs = 100)

# Generate Elbow Plots and Calculate PCs to use
Elbow <- ElbowPlot(samp, ndims = 100, reduction = "pca")
pdf(paste0(params.sampleName,"_ElbowPlot.pdf"), width = 20, height = 15)

if (params.PCs == "auto"){
    pcCount <- 1
    while(Elbow$data$stdev[pcCount]-Elbow$data$stdev[pcCount+1]>0.01 | Elbow$data$stdev[pcCount+1]-Elbow$data$stdev[pcCount+2]>0.01){
        pcCount <- pcCount + 1
    }
    params.pcMax <- pcCount
}else{
    params.pcMax <- params.PCs
}

# Find Neighbors 

samp <- FindNeighbors(samp, dims = 1:params.pcMax, reduction = "pca")

# Find Clusters
samp <- FindClusters(samp)

# Build UMAP
samp <- RunUMAP(samp, dims = 1:params.pcMax, reduction = "pca")

# Get the number of recovered cells
params.CellsRecovered <- length(readLines((paste0(params.DataDir,"/",list.files(path = params.DataDir, pattern = "barcodes")))))

#####################
# Identify Doublets #
#####################
#    (Percent Doublet is calculated based on a 0.8% increase per 1000 cells recovered as per the 10X website.)
sweep.res.list_meta <- paramSweep(samp, PCs = 1:params.pcMax, sct = F )
sweep.stats_meta <- summarizeSweep(sweep.res.list_meta, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_meta)
pK <- bcmvn_meta %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
pK <- as.numeric(as.character(pK[[1]]))
annotations <- samp$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(((params.CellsRecovered/1000)*0.008)*nrow(samp@meta.data))
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))
samp <- doubletFinder(samp,
                      PCs = 1:params.pcMax,
                      pN = 0.25,
                      pK = pK,
                      nExp = nExp_poi.adj,
                      reuse.pANN = F,
                      sct = F)


# Visualize and Count Doublets
doubletmeta <- as.character(colnames(samp[[grep("DF",(colnames(samp[[]])))]]))
pdf(paste0(params.sampleName,"VisualizeDoublets.pdf"), width = 20, height = 15)
DimPlot(samp, group.by = doubletmeta)+ ggtitle(params.sampleName)
dev.off()
DoubletCount <- table(samp[[doubletmeta]])

########################
# Subset Seurat Object #
########################
samp <- subset(samp, subset = !!as.name(doubletmeta) == "Singlet"  )

######################
# Save Seuart Object #
######################
SaveSeuratRds(samp, file = paste0(params.sampleName, "_DoubletsRemoved.rds"))

sink(paste0(params.sampleName,"validation.log"))

print(DoubletCount)
samp

sink()
