### GO ANLAYSIS 

library(org.Hs.eg.db)
library(GO.db)


go_term <- "GO:0031012"

genes <- AnnotationDbi::select(org.Hs.eg.db, 
                               keys = go_term,
                               columns = c("SYMBOL", "GENENAME"),
                               keytype = "GO")

print(genes)

write.csv(genes, file = "GOInfo/ECM_GO_0031012.csv", quote = F)

go_term <- "GO:0001837"

genes <- AnnotationDbi::select(org.Hs.eg.db, 
                               keys = go_term,
                               columns = c("SYMBOL", "GENENAME"),
                               keytype = "GO")

print(genes)

write.csv(genes, file = "GOInfo/EMT_GO_0001837.csv", quote = F)

go_term <- "GO:0016477"

genes <- AnnotationDbi::select(org.Hs.eg.db, 
                               keys = go_term,
                               columns = c("SYMBOL", "GENENAME"),
                               keytype = "GO")

print(genes)

write.csv(genes, file = "GOInfo/CM_GO_0016477.csv", quote = F)

### AFTER EXPORT, searched human genes in Axo Feature file. And reupload. 

library(Seurat)
library(viridis)
library(patchwork)

ECMgenelist <- read.csv("GOInfo/ECM_AxoGenes.tsv", sep = "\t", header = F, col.names = c("ID","SYMBOL","GE"))$SYMBOL
EMTgenelist <- read.csv("GOInfo/EMT_AxoGenes.tsv", sep = "\t", header = F, col.names = c("ID","SYMBOL","GE"))$SYMBOL
CMgenelist <- read.csv("GOInfo/CM_AxoGenes.tsv", sep = "\t", header = F, col.names = c("ID","SYMBOL","GE"))$SYMBOL

featuresSO <- as.data.frame(rownames(MergedSO))

myECMgenelist <- ECMgenelist[ECMgenelist %in% featuresSO$`rownames(MergedSO)`]
myEMTgenelist <- EMTgenelist[EMTgenelist %in% featuresSO$`rownames(MergedSO)`]
myCMgenelist <- CMgenelist[CMgenelist %in% featuresSO$`rownames(MergedSO)`]

MergedSO <- AddModuleScore(MergedSO, features = list(myECMgenelist), name = "ECM_Score")
MergedSO <- AddModuleScore(MergedSO, features = list(myEMTgenelist), name = "EMT_Score")
MergedSO <- AddModuleScore(MergedSO, features = list(myCMgenelist), name = "CM_Score")


p1 <- FeaturePlot(MergedSO, features = "ECM_Score1",  order = F, pt.size = 1) + scale_color_viridis(limits =c(min(MergedSO$ECM_Score1),max(MergedSO$ECM_Score1)), direction = -1)
p2 <- FeaturePlot(MergedSO, features = "EMT_Score1", order = F, pt.size = 1) + scale_color_viridis(limits =c(min(MergedSO$EMT_Score1),max(MergedSO$EMT_Score1)), direction = -1)
p3 <- FeaturePlot(MergedSO, features = "CM_Score1", order = F, pt.size = 1) + scale_color_viridis(limits =c(min(MergedSO$CM_Score1),max(MergedSO$CM_Score1)), direction = -1)

p4 <- VlnPlot(MergedSO, features = "ECM_Score1", group.by = "CellIdent")
p5 <- VlnPlot(MergedSO, features = "EMT_Score1", group.by = "CellIdent")
p6 <- VlnPlot(MergedSO, features = "CM_Score1", group.by = "CellIdent")

p7 <- DimPlot(MergedSO, group.by = "CellIdent", shuffle = T, pt.size = 1, label = T)  
p8 <- DimPlot(MergedSO, group.by = "TimePoint", shuffle = T, pt.size = 1, label = F)  
pdf("figures/GO_plots2.pdf",width = 20, height = 25)
p1 + p4 + p2 + p5 + p3 + p6 + p8 + p7 + plot_layout(ncol = 2, nrow = 4)
dev.off()
