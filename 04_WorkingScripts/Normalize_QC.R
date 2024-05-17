#!/usr/local/bin/Rscript


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

# Mitochondrial Genes
if (toupper(args[1]) == "NULL"){
  params.mito_genes <- "AUTO"
} else {
  params.mito_genes <- read.csv(args[1], sep = "\t", header = F)
}

# nFeature subset quantiles
params.nfeature_lower <- args[2]
params.nfeature_upper <- args[3]

# nCount subset quantiles
params.ncount_lower <- args[4]
params.ncount_upper <- args[5]

# Seurat rds file
params.RDS <- args[6]

# Sample Name
params.sample_name <- args[7]

#############
# LOAD .RDS #
#############

NameSO <- params.sample_name
assign(NameSO, LoadSeuratRds(params.RDS))

####################
# CALCULATE MT PCT #
####################

if (length(params.mito_genes) == 1 & params.mito_genes != 'AUTO'){
  params.regexs <- c(params.mito_genes[1])
}else if (params.mito_genes == "AUTO"){
  params.regexs <- c('MT-', 'mt-', 'Mt-')
} else {
  params.regexs <- c(0)
}

calcMT <- function(SO, regex) {
  if (regex[[1]] == 0){
    mito.contig.sampSpec <- intersect(params.mito_genes , rownames(SO))
    SO[["percent.mt"]] <- PercentageFeatureSet(SO, features = mito.contig.sampSpec)
  }else {
    max <- -1
    for (pat in regex){
      mito.contig <- row.names(SO)[grepl(pattern = paste0('^', pat), x = row.names(SO))]
      if (length(mito.contig) > max){
        max <- length(mito.contig)
        mito.contig.sampSpec <- mito.contig
        winner <- pat
      }
    }
    SO[["percent.mt"]] <- PercentageFeatureSet(SO, features = mito.contig.sampSpec)
  }
  return(SO)
}


assign(NameSO, calcMT(get(NameSO), params.regexs))

if (sum(get(NameSO)@meta.data$percent.mt) == 0){
  cat("ERROR: Mitochondrial genes do not match features.tsv file. \nPlease ensure gene nomenclature is matching.")
  quit()
}

########################






