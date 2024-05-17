# Takes the single sample seurat objects through a 
# process of finding doublets (Using Doublet Finder)
# and removing the predicted doublets based on the 
# expected number of cells for the sample. 

# command: Rscript DoubletFinder.R <SeuratObject>.rds <metaid> <params.vars2Regress> <data_dir>

## read in the parameters 
args <- commandArgs(trailingOnly = TRUE)

# read in .rds  
samp <- readRDS(args[1])

# create a list of all genes
all.genes <- rownames(samp)

# Scale Data
samp <- ScaleData(samp, features = all.genes ,vars.to.regress = args[3])


#### Run PCA on individual samples
samp <- RunPCA(samp, features = all.genes ,npcs = 100)

#### Generate Elbow Plots and Calculate PCs to use
Elbow <- ElbowPlot(samp, ndims = 100, reduction = "pca")
pdf(paste0(args[2],"_ElbowPlot.pdf"), width = 20, height = 15)
pcCount <- 1
while(Elbow$data$stdev[pcCount]-Elbow$data$stdev[pcCount+1]>0.01 | Elbow$data$stdev[pcCount+1]-Elbow$data$stdev[pcCount+2]>0.01){
    pcCount <- pcCount + 1
}
param.pcMax <- pcCount

##### STOPPED HERE
#### Find Neighbors

```{r}
count <- 1
for (samp in SeuratObjList){
    print(paste0("Running Find Neighbors on Sample: ", samp))
    assign(samp,FindNeighbors(get(samp), dims = 1:param.pcMax[count], reduction = "pca"))
    count <- count + 1
}
```

#### Find Clusters

```{r}
for (samp in SeuratObjList){
    print(paste0("Running Find Clusters on Sample: ", samp))
    assign(samp,FindClusters(get(samp)))
}
```

#### Build UMAP

```{r}
count <- 1
for (samp in SeuratObjList){
    print(paste0("Building UMAP on Sample: ", samp))
    assign(samp,RunUMAP(get(samp), dims = 1:param.pcMax[count], reduction = "pca"))
    count <- count + 1
}
```

### Identify Doublets

Percent Doublet is calculated based on a 0.8% increase per 1000 cells recovered as per the 10X website.

```{r}
count <- 1 
for (samp in SeuratObjList){
    print(paste0("Finding Doublets on Sample: ", samp))
    print(paste0("    Proportion of Doublet Cells: ",(params.CellsRecovered[count]/1000)*0.008))

    sweep.res.list_meta <- paramSweep(get(samp), PCs = 1:param.pcMax[count], sct = F )
    sweep.stats_meta <- summarizeSweep(sweep.res.list_meta, GT = FALSE)
    bcmvn_meta <- find.pK(sweep.stats_meta)

    print(ggplot(bcmvn_meta, aes(pK, BCmetric, group = 1))+ geom_point() + geom_line()+ ggtitle(samp))

    pK <- bcmvn_meta %>%
        filter(BCmetric == max(BCmetric)) %>%
        select(pK)

    pK <- as.numeric(as.character(pK[[1]]))
    
    print(paste0("    pK set at : ",pK))

    annotations <- get(samp)$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(((params.CellsRecovered[count]/1000)*0.008)*nrow(get(samp)@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

    assign(paste0(samp), doubletFinder(get(samp),
                                        PCs = 1:param.pcMax[count],
                                        pN = 0.25,
                                        pK = pK,
                                        nExp = nExp_poi.adj,
                                        reuse.pANN = F,
                                        sct = F))
    count <- count + 1
}
```

#### Visualize Doublets

```{r}
for (samp in SeuratObjList){
    print(paste0("Visualizing Doublets on Sample: ", samp))
    doubletmeta <- as.character(colnames(get(samp)[[grep("DF",(colnames(get(samp)[[]])))]]))
    print(DimPlot(get(samp), group.by = doubletmeta)+ ggtitle(samp)) 
}
```

#### Count Doublets

```{r}
for (samp in SeuratObjList){
    print(paste0("Counting Doublets and Singlets on Sample: ", samp))
    doubletmeta <- as.character(colnames(get(samp)[[grep("DF",(colnames(get(samp)[[]])))]]))
    print(table(get(samp)[[doubletmeta]]))
}
```

#### Subset For Singlets

```{r}
for (samp in SeuratObjList){
    print(paste0("Subsetting for Singlets on Sample: ", samp))
    doubletmeta <- as.character(colnames(get(samp)[[grep("DF",(colnames(get(samp)[[]])))]]))
    assign(samp,subset(get(samp), subset = !!as.name(doubletmeta) == "Singlet"  ))
}
```

### Merge Seurat Objects
