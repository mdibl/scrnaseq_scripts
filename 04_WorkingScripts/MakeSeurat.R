#!/usr/local/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                     Title: MakeSeurat.R                                    ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes the raw count data (features, barcodes, matrix) and generates a Seurat           ═╣
# ╠═     Object. Options to subset features are avalible.                                       ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝
# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(stringr)
library(Matrix)
library(Seurat)
library(SeuratObject)

# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
# Read in trailing arguments
args <- commandArgs(trailingOnly = TRUE)

# File path to features, barcodes, mtx directory
params.data_directory <- args[1]

# Which column to use for Seurat Obj. instantiation
params.gene_identifier <- args[2]

# File path for target gene list
params.genes_2_rm <- args[3]

# Sample Name
params.sample_name <- args[4]

# Project Name (analysis group)
params.project_name <- gsub("\\]", "", gsub("\\[", "", gsub(", ", "-", args[5])))

# Min Cells
params.min_cells <- args[6]

# Min Features
params.min_features <- args[7]

# ╔═══════════════════╗
# ╠═ Subset Features ═╣
# ╚═══════════════════╝
if (toupper(params.genes_2_rm) == "NULL"){
  if (substr(params.data_directory, nchar(params.data_directory), nchar(params.data_directory)) == "/"){
    feature_list  <- read.csv(paste0(params.data_directory, "features.tsv.gz"), sep = "\t", header = F)
    new_feat_list <- feature_list
  }else {
    feature_list  <- read.csv(paste0(params.data_directory, "/features.tsv.gz"), sep = "\t", header = F)
    new_feat_list <- feature_list
  }
}else{
  if (substr(params.data_directory, nchar(params.data_directory), nchar(params.data_directory)) == "/"){
    gene_list     <- read.csv(params.genes_2_rm, sep = "\t", header = F)
    feature_list  <- read.csv(paste0(params.data_directory, "features.tsv"), sep = "\t", header = F)
    new_feat_list <- feature_list[-c(which(feature_list %in% gene_list))]
  }else{
    gene_list     <- read.csv(params.genes_2_rm, sep = "\t", header = F)
    feature_list  <- read.csv(paste0(params.data_directory, "/features.tsv"), sep = "\t", header = F)
    new_feat_list <- feature_list[-c(which(feature_list %in% gene_list))]
  }
}

# ╔════════════════════════════════╗
# ╠═ Create 10X Object and Subset ═╣
# ╚════════════════════════════════╝
if (toupper(params.gene_identifier) == "GENE_ID"){
  params.gene_column <- 1
}else if (toupper(params.gene_identifier) == "GENE_NAME"){
  params.gene_column <- 2
}else {
  params.gene_column <- 2
}

# Create raw 10X object
Name10X <- params.sample_name
assign(Name10X, Read10X(data.dir = params.data_directory, strip.suffix = T, gene.column = params.gene_column))

# Subset 10X Object
Name10XAnnotated <- paste0(params.sample_name,"_Ann")
assign(Name10XAnnotated, get(params.sample_name)[which(rownames(get(params.sample_name)) %in% new_feat_list[,params.gene_column]),])

# ╔════════════════════════╗
# ╠═ Create Seurat Object ═╣
# ╚════════════════════════╝
NameSO <- paste0("SO_",params.sample_name)
assign(NameSO, CreateSeuratObject(counts = get(Name10XAnnotated), project = params.sample_name, min.cells = params.min_cells, min.features = params.min_features, ))

# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
SaveSeuratRds(get(NameSO), file = paste0("00_",params.sample_name, "_InitalSO.rds"))

# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
sink(paste0("00_",params.sample_name,"_InitialValidation.log"))
print("╔══════════════════════════════════════════════════════════════════════════════════════════════╗")
print("╠  MakeSeurat.R log")
print(paste0("╠  Sample: ", params.sample_name))
print("╚══════════════════════════════════════════════════════════════════════════════════════════════╝")
print(paste0("Gene Identifier: ", params.gene_identifier))
print(paste0("ProjectName: ", params.project_name))
print(paste0("Minimum Cells: ", params.min_cells))
print(paste0("Min Features: ", params.min_features))
print("")
print("Seurat Object Status:")
print(get(NameSO))
sink()

sink(paste0("00_",params.sample_name,"_InitialVersions.log"))
print("╔══════════════════════════════════════════════════════════════════════════════════════════════╗")
print("╠  MakeSeurat.R Versions")
print(paste0("╠  Sample: ", params.sample_name))
print("╚══════════════════════════════════════════════════════════════════════════════════════════════╝")
sessionInfo()
sink()
