# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                       Title: DestinyLaunch.R                                     ═╣
# ╠═                                    Updated: Sept 30 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Run Destiny Trajectory Analysis                                                        ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝

# ╔══════════════════════╗
# ╠═ Setting Parameters ═╣
# ╚══════════════════════╝
message("Setting Parameters")

# Project Name for this analysis group
params.ProjectName <- "NMT_FCT"

# Seurat Object
params.seuratObject <- readRDS(".")

params.pcMax <- "NULL"

# Resolutions
params.npcs <- c(4, 8, 12, 15, 20, 25)

# Distance Metrics
params.distancemetrics <- c('euclidean', 'cosine', 'rankcor')

source("/Users/rseaman/Desktop/016_SingleCell/99_SCpipelineScripts/scrnaseq_scripts/01_FollowOnAnalysis/Destiny.R", echo = T)
