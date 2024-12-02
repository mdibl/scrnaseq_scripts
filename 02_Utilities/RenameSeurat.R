## Rename Features

library(dplyr)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(Seurat.utils)
library(stringr)
library(patchwork)

MergedSO

# Read in new feature file
features <- read.csv("./RenameFiles/features.tsv", header = F,sep = "\t")

# Grab existing features in seurat object
inFeatures <- as.data.frame(MergedSO@assays$RNA@layers$counts@Dimnames[[1]])

# Caputre the order of the features in seurat object
inFeatures$ordering <- row.names(inFeatures)

# Zero pad the ordering (making sure enough padding) (requires stringr)
inFeatures$ordering <- str_pad(inFeatures$order, 5, pad = "0" )

# merge in new names
outFeatures <- merge(features,inFeatures,by.x = "V1",by.y = "MergedSO@assays$RNA@layers$counts@Dimnames[[1]]", all = F )

# re-order the merged table with original seurat object order
outFeatures <- outFeatures[order(outFeatures$ordering),]


# create new seurat object with new names (requires Seurat.utils package )
MergedSORenamed <- RenameGenesSeurat(MergedSO,newnames = outFeatures$V2)
?RenameGenesSeurat

# Test Rename
one <- FeaturePlot(MergedSO, features = "AMEX60DD042090", pt.size = 1, order = T)
one
two <- FeaturePlot(MergedSORenamed, features = "MPEG1::AMEX60DD007219::MSTRG.7695", pt.size = 1, order = T)
one + two
