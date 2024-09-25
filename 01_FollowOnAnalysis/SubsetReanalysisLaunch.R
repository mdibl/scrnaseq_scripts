# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                Title: SubsetReanalysisLaunch.R                             ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes a full Seurat object, subsets the object for a specific set of clusters and      ═╣
# ╠═     and reruns Scaling, PCA, Find Neighbors, Clusters, Markers, Correlation, Dimensional   ═╣
# ╠═     Reduction, and Loupe File Creation.                                                    ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝

# ╔══════════════════════╗
# ╠═ Setting Parameters ═╣
# ╚══════════════════════╝
message("Setting Parameters")

# Project Name for this analysis group
params.ProjectName <- "NMT_FCT"

# Variable To Regress For Scaling
params.VarsToRegress <- c("nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score")


# Variable Features (VF) or ALL (default VF)
params.scaleFeatures <- "VF"

# Scale Options ( SD or default SCT)
params.scaleMethod <- "SCT"

params.seuratObject <- readRDS("./RDS/NMT_LabelTransfer.rds")

params.subsetMetadataCol <- "NewCellIdent"

params.subsetIdents <- c("Neo FCT 1","Neo FCT 2","Neo FCT 3","Meta FCT 1","Meta FCT 2")

params.pcMax <- "NULL"

# Resolutions
params.Resolutions <- c(0.05,0.1,0.3,0.5,0.7,0.9,1.2,1.5)

# Integration Method: Options( CCA, RPCA, Harmony, FastMNN, NULL) where NULL is to not run
params.IntegrationMethod <- "NULL"

params.MakeLoupe <- "TRUE"

# Loupe EULA 
params.10xEULA <- "AGREE"

source("/Users/rseaman/Desktop/016_SingleCell/99_SCpipelineScripts/scrnaseq_scripts/01_FollowOnAnalysis/SubsetReanalysis.R", echo = T)
