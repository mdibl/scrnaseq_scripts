#!/usr/local/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                    Title: Normalize_QC.R                                   ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes in the Seurat Object generated from MakeSeurat.R. Calculates Cell Cycle Scores,  ═╣
# ╠═     Mitochondrial Percentages to educate a first round Quality Control. Normalization is   ═╣
# ╠═     also run at this stage.                                                                ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝
# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(dplyr)
library(stringr)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(Seurat.utils)
library(ggplot2)
library(patchwork)

# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
args <- commandArgs(trailingOnly = TRUE)

# Mitochondrial Genes
params.mito_genes <- NULL
tryCatch({
    params.mito_genes <- read.csv(args[1])$MTgenes
    params.mito_genes <- params.mito_genes[!(params.mito_genes == "" | is.na(params.mito_genes))]
}, error = function(e) {
    message("Gene List File not Provided -- Using Automatic MT detection")
})
if (length(params.mito_genes) == 0 ){
    params.mito_genes <- NULL
    message("MT Column of Gene List File is empty or does not exist -- Using Automatic MT detection")
}

# NULL means auto for logical checks
params.g2m_genes <- NULL
params.s_genes <- NULL
tryCatch({
    params.g2m_genes <- read.csv(args[1])$G2Mgenes
    params.s_genes <- read.csv(args[1])$Sgenes
    params.g2m_genes <- params.g2m_genes[!(params.g2m_genes == "" | is.na(params.g2m_genes))]
    params.s_genes <- params.s_genes[!(params.s_genes == "" | is.na(params.s_genes))]
}, error = function(e) {
    message("Gene List File not Provided -- Using Default CC detection")
})

# nFeature subset quantiles
params.nfeature_lower <- args[2]
params.nfeature_upper <- 100 - as.integer(args[3])

# nCount subset quantiles
params.ncount_lower <- args[4]
params.ncount_upper <- 100 - as.integer(args[5])

# Mito pct threshold
params.mito_pct <- args[6]

# Seurat rds file
params.SeuratObject <- args[7]

# Sample Name
params.sample_name <- args[8]

# Run CC Score
params.runCCScore <- args[9]

# ╔════════════════════╗
# ╠═ CC Scoring Logic ═╣
# ╚════════════════════╝
if (length(params.g2m_genes) == 0 | length(params.s_genes) == 0 ){
    if (length(params.g2m_genes) == length(params.s_genes)) {
        message("G2M and S gene Columns are empty or do not exist -- Using AUTO")
        params.g2m_genes <- NULL
        params.s_genes <- NULL
    }else{
        message("G2M or S gene list is empty -- Please fill both lists or leave them both empty (auto)")
        quit(status = 1)    
    }
}

if (length(params.g2m_genes) == 0 & length(params.s_genes) == 0){
    params.g2m_genes <- NULL
    params.s_genes <- NULL
    if (toupper(params.runCCScore) == "TRUE"){
        message("Cell Cycle Soring -- Auto Gene Lists")
    }else if (toupper(params.runCCScore) == "FALSE"){
        message("Cell Cycle Soring -- Skipped")
    }
}else if (length(params.g2m_genes) != 0 & length(params.s_genes) != 0){
    if (toupper(params.runCCScore) == "TRUE"){
        message("Cell Cycle Soring -- Manual Gene Lists")
    }else if (toupper(params.runCCScore) == "FALSE"){
        message("Cell Cycle Scores are not being regressed but gene list is provided")
        quit(status = 1)    
    }
}else {
    message("G2M or S gene list is empty -- Please fill both lists or leave them both empty")
    quit(status = 1)    
}


# ╔════════════════════╗
# ╠═ Load Seurat .rds ═╣
# ╚════════════════════╝
NameSO <- params.sample_name
assign(NameSO, readRDS(params.SeuratObject))

# ╔════════════════════╗
# ╠═ Calculate Mito % ═╣
# ╚════════════════════╝
if ( length(params.mito_genes) == 0){
  params.regexs <- c('MT-', 'mt-', 'Mt-')
}else {
  params.regexs <- c(params.mito_genes)
} 

calcMT <- function(SO, regex) {
    params.mito_genes
  if (length(params.mito_genes) != 0){
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
  cat("ERROR: Mitochondrial genes do not match features.tsv file or had 0 expression recorded. \nPlease ensure gene nomenclature is matching.")
  quit(status = 1)
}

# ╔═══════════════════════════════════════════╗
# ╠═ Normalize Data & Find Variable Features ═╣
# ╚═══════════════════════════════════════════╝
assign(NameSO, NormalizeData(get(NameSO)))
assign(NameSO, FindVariableFeatures(get(NameSO)))

# ╔═════════════════════════════╗
# ╠═ Generate Pre-Filter Plots ═╣
# ╚═════════════════════════════╝
plot1 <- FeatureScatter(get(NameSO), feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle = T)
plot2 <- FeatureScatter(get(NameSO), feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle = T)
vplot1 <- VlnPlot(get(NameSO), features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1)

pdf(file = paste0(params.sample_name,'_preFilterQC.pdf'), title = paste0(params.sample_name,' PreFilter QC'), width = 11, height = 8.5)
print(plot1 + plot2)
print(vplot1)
dev.off()

# ╔═══════════════════════╗
# ╠═ Subset via Metadata ═╣
# ╚═══════════════════════╝
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

# ╔══════════════════════════════╗
# ╠═ Generate Post-Filter Plots ═╣
# ╚══════════════════════════════╝
plot1 <- FeatureScatter(get(NameSO), feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle = T)
plot2 <- FeatureScatter(get(NameSO), feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle = T)
vplot1 <- VlnPlot(get(NameSO), features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1)

pdf(file = paste0(params.sample_name,'_postFilterQC.pdf'), title = paste0(params.sample_name,' PostFilter QC'), width = 11, height = 8.5)
print(plot1 + plot2)
print(vplot1)
dev.off()

# ╔════════════════════════════════════════╗
# ╠═ Cell Cycle Scoring & Rename Features ═╣
# ╚════════════════════════════════════════╝
new_features <- rownames(get(NameSO))

for (i in 1:length(new_features)) {
  gene <- new_features[i]
  gl <- unlist(strsplit(gene, ""))
  if ("_" %in% gl) {
    gene <- gsub("_", "-", gene)
  }
  new_features[i] <- gene
}

if (toupper(params.runCCScore) == "TRUE") {
    SO.rn <- RenameGenesSeurat(get(NameSO),newnames = new_features, assay = "RNA")
    assign(NameSO, SO.rn)

    if (length(params.g2m_genes) == 0 & length(params.s_genes) == 0){
        g2m <- lapply(cc.genes.updated.2019$g2m.genes, toupper)
        s <- lapply(cc.genes.updated.2019$s.genes, toupper)
    }else {
        g2m <- lapply(gsub("_","-",params.g2m_genes), toupper)
        s <- lapply(gsub("_","-",params.s_genes), toupper)
    }

    upper_gns <- lapply(new_features, toupper)

    g2m_inSO <- new_features[upper_gns %in% g2m]
    s_inSO <- new_features[upper_gns %in% s]

    assign(NameSO, CellCycleScoring(get(NameSO), s.features = s_inSO, g2m.features = g2m_inSO, set.ident = TRUE))
}

# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
SaveSeuratRds(get(NameSO), file = paste0("01_",params.sample_name, "_NormQCSO.rds"))

# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
sink(paste0("01_",params.sample_name,"_NormQCValidation.log"))
print("╔══════════════════════════════════════════════════════════════════════════════════════════════╗")
print("╠  Normalize_QC.R log")
print(paste0("╠  Sample: ", params.sample_name))
print("╚══════════════════════════════════════════════════════════════════════════════════════════════╝")
cat(paste0("\nMin nCount: ", minNCount))
cat(paste0("\nMin nFeature: ", minNFeature))
cat(paste0("\nMax nCount: ", maxNCount))
cat(paste0("\nMax nFeature: ", maxNFeature))
cat("\n")
if (toupper(params.runCCScore) == "TRUE") {
cat(paste0("\nPct G2M: " ,(length(which(get(NameSO)@meta.data$Phase == "G2M"))/length(colnames(get(NameSO))))))
cat(paste0("\nPct S: " ,(length(which(get(NameSO)@meta.data$Phase == "S"))/length(colnames(get(NameSO))))))
cat(paste0("\nPct G1: " ,(length(which(get(NameSO)@meta.data$Phase == "G1"))/length(colnames(get(NameSO))))))
}else{
    cat("\nCell Cycle Scoring was skipped")
}
cat("\n")
print("Seurat Object Status:")
print(get(NameSO))
sink()


sink(paste0("01_",params.sample_name,"_NormQCVersions.log"))
print("╔══════════════════════════════════════════════════════════════════════════════════════════════╗")
print("╠  Normalize_QC.R Versions")
print(paste0("╠  Sample: ", params.sample_name))
print("╚══════════════════════════════════════════════════════════════════════════════════════════════╝")
sessionInfo()
sink()
