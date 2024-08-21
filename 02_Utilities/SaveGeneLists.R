library(gprofiler2)
library(Seurat)

#translate <- read.table("/Users/rseaman/Desktop/000-misc/scratch/ZFgeneLists/features.tsv")
#rownames(translate) <- translate$V1

param.g2mGeneList       <- gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name # List of g2m Genes (currently set for mouse)
#param.g2mGeneList <- rownames(translate)[translate$V2 %in% param.g2mGeneList]

param.sGeneList         <- gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name # List of s Genes (currently set for mouse)
#param.sGeneList <- rownames(translate)[translate$V2 %in% param.sGeneList]

#write.csv(param.g2mGeneList, "/Users/rseaman/Desktop/000-misc/scratch/ZFgeneLists/G2Mens.tsv", row.names = F, quote = F)


#write.csv(param.sGeneList, "/Users/rseaman/Desktop/000-misc/scratch/ZFgeneLists/Sens.tsv", row.names = F, quote = F)

write.csv(param.g2mGeneList, "/Users/rseaman/Desktop/016_SingleCell/99_SCpipelineScripts/data/mouse/MM_G2M.csv", row.names = F, quote = F)


write.csv(param.sGeneList, "/Users/rseaman/Desktop/016_SingleCell/99_SCpipelineScripts/data/mouse/MM_S.csv", row.names = F, quote = F)

cc.genes.updated.2019$g2m.genes


?gorth()

cc.genes.updated.2019$g2m.genes
param.g2mGeneList
