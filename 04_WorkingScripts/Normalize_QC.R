# title: "SeuratV5_Normalize_QC.Rmd"
# author: "Ryan Seaman"
# date: "02/06/2024"

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

# File path to features, barcodes, mtx
params.data_directory <- args[1]

# Which column to use for Seurat Obj. instantiation
params.gene_identifier <- args[2]

# File path for target gene list
params.genes_2_rm <- args[3]

# Sample Name
params.sample_name <- args[4]

# Project Name (analysis group)
params.project_name <- args[5]

# Min Cells
params.min_cells <- args[6]

# Min Features
params.min_features <- args[7]

# Mitochondrial Genes
params.

######################
# FEATURE SUBSETTING #
######################

if (params.genes_2_use == null){
  feature_list  <- read.csv(paste0(params.data_directory, "features.tsv"), sep = "\t", header = F)
  new_feat_list <- feature_list
}else{
  gene_list     <- read.csv(params.genes_2_rm, sep = "\t", header = F)
  feature_list  <- read.csv(paste0(params.data_directory, "features.tsv"), sep = "\t", header = F)
  new_feat_list <- feature_list[-c(which(feature_list %in% gene_list))]
}


################################
# CREATE 10X OBJECT AND SUBSET #
################################

if (toupper(params.gene_identifier) == "GENE_ID"){
  params.gene_column <- 1
}else if (toupper(params.gene_identifier) == "GENE_NAME"){
  params.gene_column <- 2
}else {
  params.gene_column <- 2
}

# Create raw 10x object
Name10X <- params.sample_name
assign(Name10X, Read10X(data.dir = params.data_directory, strip.suffix = T, gene.column = params.gene_column))

# Subset 10x Object
Name10XAnnotated <- paste0(params.sample_name,"_Ann")
assign(Name10XAnnotated, get(params.sample_name)[which(rownames(get(params.sample_name)) %in% rownames(new_feat_list)),])

########################
# CREATE SEURAT OBJECT #
########################

NameSO <- paste0("SO_",params.sample_name)
assign(NameSO, CreateSeuratObject(counts = get(Name10XAnnotated), project = params.project_name, min.cells = params.min_cells, min.features = params.min_features))

####################
# CALCULATE MT PCT #
####################

calcMT <- function(SO) {
  mito.contig.sampSpec <- intersect(mito.contig$V2 , rownames(SO))
  SO[["percent.mt"]] <- PercentageFeatureSet(SO, features = mito.contig.sampSpec)
  SO
}

assign(NameSO, calcMT(get(NameSO)))

if (sum(get(NameSO)@meta.data$percent.mt) == 0){
  cat("ERROR: Mitochondrial genes do not match features.tsv file. \nPlease ensure gene nomenclature is matching.")
  quit()
}

########################






