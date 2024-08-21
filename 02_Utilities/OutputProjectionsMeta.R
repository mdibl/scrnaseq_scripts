

UI_UMAP <- as.data.frame(MergedSO@reductions$umap.unintegrated@cell.embeddings)
IH_UMAP <- as.data.frame(MergedSO@reductions$umap.Harmony@cell.embeddings)
UI_TSNE <- as.data.frame(MergedSO@reductions$tsne.unintegrated@cell.embeddings)
IH_TSNE <- as.data.frame(MergedSO@reductions$tsne.Harmony@cell.embeddings)
Meta <- as.data.frame(MergedSO@meta.data)

merge()

?cbind


ggplot(myT, aes(x = umapunintegrated_1, y = umapunintegrated_2)) +
    geom_point()
    
  
DimPlot(MergedSO) 
library(ggplot2)   
library(Seurat)
library(dplyr)
?try

?merge

MultiQCTable <- bind_cols(UI_UMAP,IH_UMAP,UI_TSNE,IH_TSNE,Meta)

?bind_cols

p1 <- ggplot(MultiQCTable, aes(x = umapunintegrated_1, y = umapunintegrated_2, colour = orig.ident)) +
    geom_point()


p2 <- DimPlot(MergedSO, group.by = "orig.ident") 

?ggplot
library(patchwork)


p1 + p2

test <- "l;kasjdfl;ksjdfl;askjfas"


print(test)
cat(test)
message(test)
