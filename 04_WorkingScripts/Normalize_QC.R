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
  params.mito_genes <- as.data.frame(read.csv(args[1], sep = "\t", header = F))
}

# nFeature subset quantiles
params.nfeature_lower <- args[2]
params.nfeature_upper <- 100 - as.integer(args[3])



# nCount subset quantiles
params.ncount_lower <- args[4]
params.ncount_upper <- 100 - as.integer(args[5])

# Mito pct threshold
params.mito_pct <- args[6]

# Seurat rds file
params.RDS <- args[7]

# Sample Name
params.sample_name <- args[8]

#############
# LOAD .RDS #
#############

NameSO <- params.sample_name
assign(NameSO, LoadSeuratRds(params.RDS))

####################
# CALCULATE MT PCT #
####################

if (length(params.mito_genes[,1]) == 1 && params.mito_genes != 'AUTO'){
  params.regexs <- c(params.mito_genes[1])
}else if ("AUTO" %in% params.mito_genes){
  params.regexs <- c('MT-', 'mt-', 'Mt-')
} else {
  params.regexs <- c(0)
}

calcMT <- function(SO, regex) {
  if (regex[[1]] == 0){
    mito.contig.sampSpec <- intersect(params.mito_genes[,1] , rownames(SO))
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
  cat("ERROR: Mitochondrial genes do not match features.tsv file or had 0 expression recorded. \nPlease ensure gene nomenclature is matching.")
  quit(status = 1)
}

############################################
# NORMALIZE DATA  & FIND VARIABLE FEATURES #
############################################

assign(NameSO, NormalizeData(get(NameSO)))
assign(NameSO, FindVariableFeatures(get(NameSO)))

####################
# PRE-FILTER PLOTS #
####################

plot1 <- FeatureScatter(get(NameSO), feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle = T)
plot2 <- FeatureScatter(get(NameSO), feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle = T)
vplot1 <- VlnPlot(get(NameSO), features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1)

pdf(file = paste0(params.sample_name,'_preFilterQC.pdf'), title = paste0(params.sample_name,' PreFilter QC'), width = 11, height = 8.5)
print(plot1 + plot2)
print(vplot1)
dev.off()

############################
# SUBSETTING VIA META DATA #
############################

minNCount         <- quantile(get(NameSO)$nCount_RNA, probs = seq(0,1,0.05))[[paste0(params.ncount_lower, "%")]]
minNFeature       <- quantile(get(NameSO)$nFeature_RNA, probs = seq(0,1,0.05))[[paste0(params.nfeature_lower, "%")]]
maxNCount         <- quantile(get(NameSO)$nCount_RNA, probs = seq(0,1,0.05))[[paste0(params.ncount_upper, "%")]]
maxNFeature       <- quantile(get(NameSO)$nFeature_RNA, probs = seq(0,1,0.05))[[paste0(params.nfeature_upper, "%")]]
maxMitoPct        <- as.integer(params.mito_pct)

print(minNFeature)
print(minNCount)
print((get(NameSO)@meta.data$nCount)[(get(NameSO)@meta.data$nCount) < 3875])
print(maxMitoPct)
assign(NameSO,subset(get(NameSO), subset = nFeature_RNA > minNFeature & nCount_RNA > minNCount & percent.mt < maxMitoPct & nFeature_RNA < maxNFeature & nCount_RNA < maxNCount))

####################
# POST-FILTER PLOTS #
####################

plot1 <- FeatureScatter(get(NameSO), feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle = T)
plot2 <- FeatureScatter(get(NameSO), feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle = T)
vplot1 <- VlnPlot(get(NameSO), features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1)

pdf(file = paste0(params.sample_name,'_postFilterQC.pdf'), title = paste0(params.sample_name,' PostFilter QC'), width = 11, height = 8.5)
print(plot1 + plot2)
print(vplot1)
dev.off()

#########################################
# CELL CYCLE SCORING & FEATURE RENAMING #
#########################################

new_features <- rownames(get(NameSO))

for (i in 1:length(new_features)) {
  
  gene <- new_features[i]
  gl <- unlist(strsplit(gene, ""))
  if ("_" %in% gl) {
    gene <- gsub("_", "-", gene)
  }
  new_features[i] <- gene
}

SO.rn <- RenameGenesSeurat(get(NameSO),newnames = new_features, assay = "RNA")
assign(NameSO, SO.rn)

g2m <- lapply(cc.genes.updated.2019$g2m.genes, toupper)
s <- lapply(cc.genes.updated.2019$s.genes, toupper)

upper_gns <- lapply(new_features, toupper)

g2m_inSO <- new_features[upper_gns %in% g2m]
s_inSO <- new_features[upper_gns %in% s]
assign(NameSO, CellCycleScoring(get(NameSO), s.features = s_inSO, g2m.features = g2m_inSO, set.ident = TRUE))



###############
# SAVE OUTPUT #
###############

SaveSeuratRds(get(NameSO), file = paste0(params.sample_name, "_QC.rds"))

sink(paste0(params.sample_name,".validation.log"))

print(get(NameSO))
cat("\n")
cat(paste0("\nMin nCount: ", minNCount))
cat(paste0("\nMin nFeature: ", minNFeature))
cat(paste0("\nMax nCount: ", maxNCount))
cat(paste0("\nMax nFeature: ", maxNFeature))
cat("\n")
cat(paste0("\nPct G2M: " ,(length(which(get(NameSO)@meta.data$Phase == "G2M"))/length(colnames(get(NameSO))))))
cat(paste0("\nPct S: " ,(length(which(get(NameSO)@meta.data$Phase == "S"))/length(colnames(get(NameSO))))))
cat(paste0("\nPct G1: " ,(length(which(get(NameSO)@meta.data$Phase == "G1"))/length(colnames(get(NameSO))))))

sink()

