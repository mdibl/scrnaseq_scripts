# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                   Title: SubsetReanalysis.R                                 ═╣
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

# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(dplyr)
library(stringr)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(ggplot2)
library(ggsankey)
library(ggcorrplot)
library(ggdendro)
library(viridis)
library(findPC)
library(presto)

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

# ╔═══════════════════╗
# ╠═ Create Log File ═╣
# ╚═══════════════════╝
sink(file = paste0(params.ProjectName,"SubsetExecution.log"), append = T)
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  Subset Execution Log\n")
cat(paste0("╠  Subset Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("╠═══════════════════════════════════  Parameter Settings  ═════════════════════════════════════╣\n")
cat(paste0("params.ProjectName: ", params.ProjectName,"\n"))
cat(paste0("params.VarsToRegress: \n"))
cat(paste0("\t",params.VarsToRegress,"\n"))
cat(paste0("params.scaleFeatures: ", params.scaleFeatures,"\n"))
cat(paste0("params.scaleMethod: ", params.scaleMethod,"\n"))
cat(paste0("params.seuratObject: \n"))
print(params.seuratObject)
cat(paste0("params.subsetMetadataCol: ", params.subsetMetadataCol,"\n"))
cat(paste0("params.subsetIdents: \n"))
cat(paste0("\t",params.subsetIdents,"\n"))
cat(paste0("params.pcMax: ", params.pcMax,"\n"))
cat(paste0("params.Resolutions: \n"))
cat(paste0("\t",params.Resolutions,"\n"))
cat(paste0("params.IntegrationMethod: ", params.IntegrationMethod,"\n"))
cat(paste0("params.MakeLoupe: ", params.MakeLoupe,"\n"))
cat(paste0("params.10xEULA: ", params.10xEULA,"\n"))
sink()

# ╔═════════════╗
# ╠═ Subset SO ═╣
# ╚═════════════╝
SubsetSO <- params.seuratObject
SubsetSO$Subset <- SubsetSO[[params.subsetMetadataCol]]
for(i in params.subsetIdents){
    SubsetSO$Subset <- sub(i,"ToBeSubset",SubsetSO$Subset)
}
sink(file = paste0(params.ProjectName,"SubsetExecution.log"), append = T)
cat("\n")
cat("╠═══════════════════════════════════════  Subsetting ══════════════════════════════════════════╣\n")
cat(paste0("Subsetting out ",sum(SubsetSO$Subset == "ToBeSubset"), " / ", nrow(SubsetSO), " cells.\n"))
sink()

message(paste0("Subsetting out ",sum(SubsetSO$Subset == "ToBeSubset"), " / ", nrow(SubsetSO), " cells."))

SubsetSO <- subset(SubsetSO, subset= Subset == "ToBeSubset")

# ╔══════════════╗
# ╠═ Cleanup SO ═╣
# ╚══════════════╝
SubsetSO@assays$SCT <- NULL
DefaultAssay(object = SubsetSO) <- "RNA"
SubsetSO$Subset <- NULL

# ╔══════════════════╗
# ╠═ Normalize Data ═╣
# ╚══════════════════╝
message("Normalizing Data")
SubsetSO <- NormalizeData(SubsetSO)

# ╔══════════════════════════════╗
# ╠═ Scale Merged Seurat Object ═╣
# ╚══════════════════════════════╝
message("Scaling Merged Seurat Object")
# create a list of all genes
all.genes <- rownames(SubsetSO)

# Scale Data and PCA
if (params.scaleMethod == "SD"){
    if (params.scaleFeatures == "ALL"){
        SubsetSO <- FindVariableFeatures(SubsetSO)
        SubsetSO <- ScaleData(SubsetSO, features = all.genes ,vars.to.regress = params.VarsToRegress)
    }else{
        SubsetSO <- FindVariableFeatures(SubsetSO)
        SubsetSO <- ScaleData(SubsetSO, vars.to.regress = params.VarsToRegress)
    }
}else{
    if (params.scaleFeatures == "ALL"){
        SubsetSO <- SCTransform(SubsetSO, vars.to.regress = params.VarsToRegress, return.only.var.genes = F )
    }else{
        SubsetSO <- SCTransform(SubsetSO, vars.to.regress = params.VarsToRegress)
    }
}

sink(file = paste0(params.ProjectName,"SubsetExecution.log"), append = T)
cat(paste0("Subset Seurat Object: \n"))
print(SubsetSO)
sink()

# ╔═══════════╗
# ╠═ Run PCA ═╣
# ╚═══════════╝
message("Running PCA")
SubsetSO <- RunPCA(SubsetSO, npcs = 100)

# ╔═══════════════╗
# ╠═ Find PC Max ═╣
# ╚═══════════════╝
message("Finding PC Max")
Elbow <- ElbowPlot(SubsetSO,  ndims = 100, reduction = "pca")
ElbowPoints <- Elbow$data

for (i in seq(0.1, 0.9, 0.1)){
    loess <- loess(stdev ~ dims,data=ElbowPoints, span = i)
    ElbowPoints[paste0("loessS", i)] <- loess$fitted
}

for (i in seq(0.1, 0.9, 0.1)){
    first_deriv <- diff(ElbowPoints[,paste0("loessS", i)])/diff(ElbowPoints$dims)
    ElbowPoints[3:100,paste0("stddev_Der2_", i)]<- (diff(first_deriv)/diff(ElbowPoints$dims))[1:length(diff(first_deriv)/diff(ElbowPoints$dims))-1]
}

columns <- colnames(ElbowPoints)
idents <- seq(0.1, 0.9, 0.1)
SSE <- c()
deriv2_var <- c()
for (item in columns){
    if (length(grep("Der2", item)) > 0){
        variance <- var(ElbowPoints[item], na.rm = T)
        deriv2_var <- append(deriv2_var, variance)
    }else if (length(grep("loess", item)) > 0 & !(length(grep("Der", item)) > 0)){
        err <- ElbowPoints[item] - ElbowPoints$stdev
        norms <- sum(((err - min(err))**2)/((max(err) - min(err))**2))
        SSE <- append(SSE, norms)
    }
}

rm_ind <- which(deriv2_var > 0.010)
if (length(rm_ind) > 0){
    SSE <- SSE[-rm_ind]
    deriv2_var <- deriv2_var[-rm_ind]
    idents <- idents[-rm_ind]
    index <- which(SSE == min(SSE))
    ident <- paste0('loessS', idents[index])
}else {
    index <- 1
    ident <- paste0('loessS', idents[index])
}

if (toupper(params.pcMax) == "NULL"){
    for (i in 1:length(ElbowPoints[[ident]])) {
        tmp <- ElbowPoints[[ident]][i]
        if (i == length(ElbowPoints[[ident]])){
            break
        }
        if (ElbowPoints[[ident]][i+1] >= tmp){
            ElbowPoints[[ident]][i+1] <- tmp - 0.00001
        }
    }
    
    pc_tbl <- findPC(sdev = ElbowPoints[[ident]], number = 100, method = "all", figure = T)
    params.pcMax <- mean(x = c(pc_tbl[1,2], pc_tbl[1,3], pc_tbl[1,4]))
    params.pcMax <- ceiling(params.pcMax)
}else{
    params.pcMax
}

x = ElbowPoints$dims
y = ElbowPoints$stdev
df <- data.frame(x,y)

pdf(paste0(params.ProjectName,"_Merged_ElbowPlot.pdf"), width = 20, height = 15)
ggplot(df, aes(x, y)) +
    geom_point() +
    geom_line(aes(x = seq(1,100), y = ElbowPoints[[ident]]), color = "green") +
    xlab("PC") +
    ylab("Std Dev") +
    geom_vline(xintercept = params.pcMax, linetype="dotted", 
               color = "red", size=1.5) +
    ggtitle(paste0("Loess Regression of Std Dev ~ PC  :  PC Chosen = ", params.pcMax))
dev.off()

sink(file = paste0(params.ProjectName,"SubsetExecution.log"), append = T)
cat("\n")
cat("╠═══════════════════════════════════════════  PCA  ════════════════════════════════════════════╣\n")
cat("Loess Summary: \n")
print(my_summary)
cat("\n")
cat("PC cutoff calculations:\n")
cat(paste0("First Derivative : ", pc_tbl[1,2], "\n"))
cat(paste0("Second Derivative : ", pc_tbl[1,3], "\n"))
cat(paste0("Preceding Residual : ", pc_tbl[1,4], "\n"))
cat("\n")
cat(paste0("PCs used: 1 - ", params.pcMax,"\n"))
sink()

# ╔═══════════════════╗
# ╠═ Run Integration ═╣
# ╚═══════════════════╝
message("Running Integration")
if(params.IntegrationMethod == "FastMNN"){
    if(params.scaleMethod == "SD"){
        SubsetSO <- IntegrateLayers(object = SubsetSO, method = paste0(params.IntegrationMethod,"Integration"), new.reduction = paste0("integrated.",params.IntegrationMethod))
    }else{
        SubsetSO <- IntegrateLayers(object = SubsetSO, method = paste0(params.IntegrationMethod,"Integration"), new.reduction = paste0("integrated.",params.IntegrationMethod), normalization.method = "SCT")
    }
}else if (params.IntegrationMethod == "NULL"){
    message("Not Integrating")
}else{
    if(params.scaleMethod == "SD"){
        SubsetSO <- IntegrateLayers(object = SubsetSO, method = paste0(params.IntegrationMethod,"Integration"), orig.reduction = "pca", new.reduction = paste0("integrated.",params.IntegrationMethod))
    }else{
        SubsetSO <- IntegrateLayers(object = SubsetSO, method = paste0(params.IntegrationMethod,"Integration"), orig.reduction = "pca", new.reduction = paste0("integrated.",params.IntegrationMethod), normalization.method = "SCT")
    }
}

sink(file = paste0(params.ProjectName,"SubsetExecution.log"), append = T)
cat("\n")
cat("╠═══════════════════════════════════════  Integration  ════════════════════════════════════════╣\n")
cat(paste0("Integrated Seurat Object: \n"))
print(SubsetSO)
sink()


# ╔══════════════════════╗
# ╠═ Create Directories ═╣
# ╚══════════════════════╝
message("Creating Marker Directories")
dir.create("markers")
dir.create("markers/unintegrated/")
if (params.IntegrationMethod != "NULL"){
    dir.create("markers/integrated/")
}

# ╔════════════════════════════════════╗
# ╠═ Find Neighbors and Find Clusters ═╣
# ╚════════════════════════════════════╝
message("Finding Neighbors and Clusters")
# Generate Resolution Names -- unintegrated
count <- 1
countMax <- length(params.Resolutions)
params.UnintegratedResolutionsNames <- c()
while (count <= countMax){
    resVal <- paste0("unintegratedSubsetRes.",params.Resolutions[count])
    params.UnintegratedResolutionsNames <- append(params.UnintegratedResolutionsNames, resVal)
    count <- count + 1
}
params.UnintegratedResolutionsNames

# Find Neighbors -- unintegrated
SubsetSO <- FindNeighbors(SubsetSO, dims = 1:params.pcMax, reduction = "pca")

# Find Clusters -- unintegrated
SubsetSO <- FindClusters(SubsetSO, resolution = params.Resolutions, cluster.name = params.UnintegratedResolutionsNames)

if (params.IntegrationMethod != "NULL"){
    # Generate Resolution Names -- integrated
    count <- 1
    countMax <- length(params.Resolutions)
    params.IntegratedResolutionsNames <- c()
    while (count <= countMax){
        resVal <- paste0(params.IntegrationMethod,"SubsetRes.",params.Resolutions[count])
        params.IntegratedResolutionsNames <- append(params.IntegratedResolutionsNames, resVal)
        count <- count + 1
    }
    params.IntegratedResolutionsNames
    # Find Neighbors -- integrated
    SubsetSO <- FindNeighbors(SubsetSO, dims = 1:params.pcMax, reduction = paste0("integrated.",params.IntegrationMethod))
    # Find Cluster -- integrated
    SubsetSO <- FindClusters(SubsetSO, resolution = params.Resolutions, cluster.name = params.IntegratedResolutionsNames)
}

# ╔════════════════════════╗
# ╠═ Generate Sankey Plot ═╣
# ╚════════════════════════╝
message("Generating Sankey Plots")
SubsetSOmetaRes <- SubsetSO@meta.data
SubsetSOmetaRes <- SubsetSOmetaRes[,grepl("^unintegratedSubsetRes",names(SubsetSOmetaRes))]
names(SubsetSOmetaRes) <- gsub("unintegratedSubsetRes","res",names(SubsetSOmetaRes))

df <- SubsetSOmetaRes%>%
    make_long(names(SubsetSOmetaRes))

df$next_node <- str_pad(df$next_node,2, pad = "0")
df$node <- str_pad(df$node,2, pad = "0")

df$reswClust <- paste0(df$x,": ", df$node)

reagg <- df%>%
    dplyr::group_by(reswClust)%>%                                        # Get counts 
    tally()

df2 <- merge(df, 
             reagg, 
             by.x = 'reswClust', 
             by.y = 'reswClust', 
             all.x = TRUE)

pl <- ggplot(df2, aes(x = x,                        
                      next_x = next_x,                                     
                      node = node,
                      next_node = next_node,        
                      fill = factor(node),
                      label = paste0(node, " = ", n))) +                 # Creates lable with node and count
    geom_sankey(flow.alpha = 0.5,                                        # Add opacity
                node.color = "black",                                    # Node Color
                show.legend = TRUE) +                                    # This determines if you want your legend to show
    geom_sankey_label(size = 3,
                      color = "black", 
                      fill = "white") +                                  # This specifies the Label format for each node 
    theme_bw() + 
    theme(legend.position = 'none') + 
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) + 
    scale_fill_viridis_d(option = "H") + 
    labs(title = paste0("Sankey: ",params.ProjectName, " Unintegrated Subset Resolutions")) + 
    labs(caption =NULL ) + 
    labs(fill = 'Nodes')

### Save Sankey Plot
pdf(paste0(params.ProjectName,"UnintegratedSankeyPlot.pdf"), width = 20, height = 15)
pl
dev.off()


if (params.IntegrationMethod != "NULL"){
    SubsetSOmetaRes <- SubsetSO@meta.data
    pattern <- paste0(params.IntegrationMethod,"SubsetRes")
    SubsetSOmetaRes <- SubsetSOmetaRes[,grepl(pattern,names(SubsetSOmetaRes))]
    names(SubsetSOmetaRes) <- gsub(pattern,"res",names(SubsetSOmetaRes))
    
    df <- SubsetSOmetaRes%>%
        make_long(names(SubsetSOmetaRes))
    
    df$next_node <- str_pad(df$next_node,2, pad = "0")
    df$node <- str_pad(df$node,2, pad = "0")
    
    df$reswClust <- paste0(df$x,": ", df$node)
    
    reagg <- df%>%
        dplyr::group_by(reswClust)%>%                                        # Get counts 
        tally()
    
    df2 <- merge(df, 
                 reagg, 
                 by.x = 'reswClust', 
                 by.y = 'reswClust', 
                 all.x = TRUE)
    
    pl <- ggplot(df2, aes(x = x,                        
                          next_x = next_x,                                     
                          node = node,
                          next_node = next_node,        
                          fill = factor(node),
                          label = paste0(node, " = ", n))) +                 # Creates lable with node and count
        geom_sankey(flow.alpha = 0.5,                                        # Add opacity
                    node.color = "black",                                    # Node Color
                    show.legend = TRUE) +                                    # This determines if you want your legend to show
        geom_sankey_label(size = 3,
                          color = "black", 
                          fill = "white") +                                  # This specifies the Label format for each node 
        theme_bw() + 
        theme(legend.position = 'none') + 
        theme(axis.title = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank()) + 
        scale_fill_viridis_d(option = "H") + 
        labs(title = paste0("Sankey: ",params.ProjectName," ", params.IntegrationMethod," Integrated Subset Resolutions")) + 
        labs(caption =NULL ) + 
        labs(fill = 'Nodes')
    
    ### Save Sankey Plot
    pdf(paste0(params.ProjectName,params.IntegrationMethod ,"SankeyPlot.pdf"), width = 20, height = 15)
    print(pl)
    dev.off()
}
# ╔═══════════════╗
# ╠═ Join Layers ═╣
# ╚═══════════════╝
message("Joining Layers")
SubsetSO[["RNA"]] <- JoinLayers(SubsetSO[["RNA"]])

# ╔════════════════════════════════════╗
# ╠═ Find Markers and Save TSV output ═╣
# ╚════════════════════════════════════╝
message("Finding Markers and Saving out Marker Files")
if (params.scaleMethod != "SD"){
    SubsetSO <- PrepSCTFindMarkers(SubsetSO, assay = "SCT", verbose = TRUE)
}

for(i in 1:length(params.Resolutions)){
    print(paste0("Finding Unintegrated Markers at ",params.Resolutions[i]," resolution..."))
    Idents(SubsetSO) <- paste0("unintegratedSubsetRes.",params.Resolutions[i])
    assign(paste0("UnintegratedMarkersRes",params.Resolutions[i]),FindAllMarkers(SubsetSO, only.pos = T, logfc.threshold = 0.3))
}

# Add in annotation info and output tsv files
for(j in 1:length(params.Resolutions)){
    print(paste0("Generating Unintegrated Markers tsv at ",params.Resolutions[j]," resolution..."))
    markertemp <- get(paste0("UnintegratedMarkersRes",params.Resolutions[j]))
    markertemp <- markertemp  %>% arrange(cluster , desc(avg_log2FC))
    markertemp$markercall <- 1
    markertemp <- markertemp %>%
        mutate(markercall = ifelse(avg_log2FC >= 1.7 & markertemp$pct.1/markertemp$pct.2 >= 1.7 & -log10(markertemp$p_val_adj) >= 27, 1, 0))
    assign(paste0("UnintegratedMarkersRes",params.Resolutions[j]),markertemp)
    write.table(get(paste0("UnintegratedMarkersRes",params.Resolutions[j])), file = paste0("markers/unintegrated/UnintegratedMarkersRes",params.Resolutions[j],".tsv"), sep = "\t", quote = F,row.names = F)
    rm(markertemp)
}

if (params.IntegrationMethod != "NULL") {
    for(i in 1:length(params.Resolutions)){
        print(paste0("Finding ", params.IntegrationMethod ," Markers at ",params.Resolutions[i]," resolution..."))
        Idents(SubsetSO) <- paste0(params.IntegrationMethod,"SubsetRes.",params.Resolutions[i])
        assign(paste0("IntegratedMarkersRes",params.Resolutions[i]),FindAllMarkers(SubsetSO, only.pos = T, logfc.threshold = 0.3))
    }
    
    # Add in annotation info and output tsv files
    for(j in 1:length(params.Resolutions)){
        print(paste0("Generating ", params.IntegrationMethod," Markers tsv at ",params.Resolutions[j]," resolution..."))
        markertemp <- get(paste0("IntegratedMarkersRes",params.Resolutions[j]))
        markertemp <- markertemp  %>% arrange(cluster , desc(avg_log2FC))
        markertemp$markercall <- 1
        markertemp <- markertemp %>%
            mutate(markercall = ifelse(avg_log2FC >= 1.7 & markertemp$pct.1/markertemp$pct.2 >= 1.7 & -log10(markertemp$p_val_adj) >= 27, 1, 0))
        assign(paste0("IntegratedMarkersRes",params.Resolutions[j]),markertemp)
        write.table(get(paste0("IntegratedMarkersRes",params.Resolutions[j])), file = paste0("markers/integrated/",params.IntegrationMethod,"MarkersRes",params.Resolutions[j],".tsv"), sep = "\t", quote = F,row.names = F)
        rm(markertemp)
    }
}

# ╔═══════════════════╗
# ╠═ Run Correlation ═╣
# ╚═══════════════════╝
message("Running Correlation, Generating Plot, and Modifying Cluster Order in Seurat Object")
VarGenes <- VariableFeatures(SubsetSO)
pdf(paste0(params.ProjectName,"UnintegratedCorPlot.pdf"), width = 20, height = 12)
for(i in 1:length(params.Resolutions)){
    AvgExp <- as.data.frame(AggregateExpression(SubsetSO, features = VarGenes, group.by = paste0("unintegratedSubsetRes.",params.Resolutions[i]))$RNA)
    colnames(AvgExp) <- gsub("^g","Clust: ", colnames(AvgExp))
    Cor <- cor(AvgExp, method = "pearson")
    Dist <- as.dist(1 - Cor)
    DistHClust <- hclust(Dist)
    Order <- DistHClust$labels[DistHClust$order]
    Cor <- Cor[Order,]
    Cor <- Cor[,Order]
    p1 <- ggcorrplot(Cor) + scale_fill_gradientn(colors = viridis(256, option = 'D' ,direction = -1)) + ggtitle(paste0("unintegratedSubsetRes.",params.Resolutions[i]))
    p2 <- ggdendrogram(DistHClust, rotate = F)
    print(p1 + p2)
    Order <- gsub("^Clust: ","", Order)
    SubsetSO@meta.data[[paste0("unintegratedSubsetRes.",params.Resolutions[i])]] <- factor(SubsetSO@meta.data[[paste0("unintegratedSubsetRes.",params.Resolutions[i])]],levels = Order)
}
dev.off()

pdf(paste0(params.ProjectName,params.IntegrationMethod ,"CorPlot.pdf"), width = 20, height = 12)
if (params.IntegrationMethod != "NULL"){
    for(i in 1:length(params.Resolutions)){
        AvgExp <- as.data.frame(AggregateExpression(SubsetSO, features = VarGenes, group.by = paste0(params.IntegrationMethod,"SubsetRes.",params.Resolutions[i]))$RNA)
        colnames(AvgExp) <- gsub("^g","Clust: ", colnames(AvgExp))
        Cor <- cor(AvgExp, method = "pearson")
        Dist <- as.dist(1 - Cor)
        DistHClust <- hclust(Dist)
        Order <- DistHClust$labels[DistHClust$order]
        Cor <- Cor[Order,]
        Cor <- Cor[,Order]
        p1 <- ggcorrplot(Cor) + scale_fill_gradientn(colors = viridis(256, option = 'D' ,direction = -1)) + ggtitle(paste0(params.IntegrationMethod,"SubsetRes.",params.Resolutions[i]))
        p2 <- ggdendrogram(DistHClust, rotate = F)
        print(p1 + p2)
        Order <- gsub("^Clust: ","", Order)
        SubsetSO@meta.data[[paste0(params.IntegrationMethod,"SubsetRes.",params.Resolutions[i])]] <- factor(SubsetSO@meta.data[[paste0(params.IntegrationMethod,"SubsetRes.",params.Resolutions[i])]],levels = Order )
    }
}
dev.off()

if(params.IntegrationMethod == "NULL"){
    file.remove(paste0(params.ProjectName,params.IntegrationMethod ,"CorPlot.pdf"))
}

sink(file = paste0(params.ProjectName,"SubsetExecution.log"), append = T)
cat("\n")
cat("╠═══════════════════════════════════════  Clustering  ═════════════════════════════════════════╣\n")
cat(paste0("Clustered Seurat Object: \n"))
print(SubsetSO)
sink()


# ╔══════════════════════╗
# ╠═ Make UMAP and TSNE ═╣
# ╚══════════════════════╝
message("Running UMAP and TSNE Reductions")
SubsetSO <- RunUMAP(SubsetSO, dims = 1:params.pcMax, reduction = "pca", reduction.name = "umap.unintegrated.Subset")
if(params.IntegrationMethod != "NULL" ){
    SubsetSO <- RunUMAP(SubsetSO, dims = 1:params.pcMax, reduction = paste0("integrated.",params.IntegrationMethod), reduction.name = paste0("umap.",params.IntegrationMethod,".Subset"))
}
SubsetSO <- RunTSNE(SubsetSO, dims = 1:params.pcMax, reduction = "pca", reduction.name = "tsne.unintegrated.Subset")
if(params.IntegrationMethod != "NULL" ){
    SubsetSO <- RunTSNE(SubsetSO, dims = 1:params.pcMax, reduction = paste0("integrated.",params.IntegrationMethod), reduction.name = paste0("tsne.",params.IntegrationMethod,".Subset"))
}

# Set Front Page Layout
Page1layout <- "
AAAAAA
AAAAAA
AAAAAA
BBCCDD
"

# ╔══════════════════════════╗
# ╠═ Plot Unintegrated UMAP ═╣
# ╚══════════════════════════╝
message("Generating Unintegrated UMAP Plots")
p1 <- DimPlot(object = SubsetSO, reduction = 'umap.unintegrated.Subset', pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = "orig.ident", shuffle = T)
p2 <- FeaturePlot(SubsetSO, reduction = "umap.unintegrated.Subset", pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "nFeature_RNA", order = T) + scale_color_viridis(limits =c(min(SubsetSO$nFeature_RNA),max(SubsetSO$nFeature_RNA)), direction = -1)
p3 <- FeaturePlot(SubsetSO, reduction = "umap.unintegrated.Subset", pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "nCount_RNA", order = T) + scale_color_viridis(limits =c(min(SubsetSO$nCount_RNA),max(SubsetSO$nCount_RNA)), direction = -1)
p4 <- FeaturePlot(SubsetSO, reduction = "umap.unintegrated.Subset", pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "percent.mt", order = T) + scale_color_viridis(limits =c(min(SubsetSO$percent.mt),max(SubsetSO$percent.mt)), direction = -1)

pdf(paste0(params.ProjectName,"UnintegratedUMAP.pdf"),width = 20, height = 15)
p1 + p2 + p3 + p4 + plot_layout(design = Page1layout)
try(expr = {DimPlot(object = SubsetSO, reduction = 'umap.unintegrated.Subset', pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = "CIscCATCH")})
for (i in params.Resolutions){
    print(DimPlot(object = SubsetSO, reduction = 'umap.unintegrated.Subset', pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = paste0("unintegratedSubsetRes.",i), shuffle = T))
}
dev.off()

# ╔══════════════════════════╗
# ╠═ Plot Unintegrated TSNE ═╣
# ╚══════════════════════════╝
message("Generating Unintegrated TSNE Plots")
p1 <- DimPlot(object = SubsetSO, reduction = 'tsne.unintegrated.Subset', pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = "orig.ident", shuffle = T)
p2 <- FeaturePlot(SubsetSO, reduction = "tsne.unintegrated.Subset", pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "nFeature_RNA", order = T) + scale_color_viridis(limits =c(min(SubsetSO$nFeature_RNA),max(SubsetSO$nFeature_RNA)), direction = -1)
p3 <- FeaturePlot(SubsetSO, reduction = "tsne.unintegrated.Subset", pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "nCount_RNA", order = T) + scale_color_viridis(limits =c(min(SubsetSO$nCount_RNA),max(SubsetSO$nCount_RNA)), direction = -1)
p4 <- FeaturePlot(SubsetSO, reduction = "tsne.unintegrated.Subset", pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "percent.mt", order = T) + scale_color_viridis(limits =c(min(SubsetSO$percent.mt),max(SubsetSO$percent.mt)), direction = -1)

pdf(paste0(params.ProjectName,"UnintegratedTSNE.pdf"),width = 20, height = 15)
p1 + p2 + p3 + p4 + plot_layout(design = Page1layout)
try(expr = {DimPlot(object = SubsetSO, reduction = 'tsne.unintegrated.Subset', pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = "CIscCATCH")})
for (i in params.Resolutions){
    print(DimPlot(object = SubsetSO, reduction = 'tsne.unintegrated.Subset', pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = paste0("unintegratedSubsetRes.",i), shuffle = T))
}
dev.off()

# ╔═════════════════════════════════╗
# ╠═ Plot Integrated UMAP (if run) ═╣
# ╚═════════════════════════════════╝
pdf(paste0(params.ProjectName,params.IntegrationMethod,"IntegratedUMAP.pdf"),width = 20, height = 15)
if (params.IntegrationMethod != "NULL"){
    message("Generating Integrated UMAP Plots")
    p1 <- DimPlot(object = SubsetSO, reduction = paste0("umap.",params.IntegrationMethod,".Subset"), pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = "orig.ident", shuffle = T)
    p2 <- FeaturePlot(SubsetSO, reduction = paste0("umap.",params.IntegrationMethod,".Subset"), pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "nFeature_RNA", order = T) + scale_color_viridis(limits =c(min(SubsetSO$nFeature_RNA),max(SubsetSO$nFeature_RNA)), direction = -1)
    p3 <- FeaturePlot(SubsetSO, reduction = paste0("umap.",params.IntegrationMethod,".Subset"), pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "nCount_RNA", order = T) + scale_color_viridis(limits =c(min(SubsetSO$nCount_RNA),max(SubsetSO$nCount_RNA)), direction = -1)
    p4 <- FeaturePlot(SubsetSO, reduction = paste0("umap.",params.IntegrationMethod,".Subset"), pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "percent.mt", order = T) + scale_color_viridis(limits =c(min(SubsetSO$percent.mt),max(SubsetSO$percent.mt)), direction = -1)
    
    print(p1 + p2 + p3 + p4 + plot_layout(design = Page1layout))
    try(expr = {print(DimPlot(object = SubsetSO, reduction = paste0("umap.",params.IntegrationMethod,".Subset"), pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = "CIscCATCH"))})
    for (i in params.Resolutions){
        print(DimPlot(object = SubsetSO, reduction = paste0("umap.",params.IntegrationMethod,".Subset"), pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = paste0(params.IntegrationMethod,"SubsetRes.",i)))
    }
}
dev.off()

if(params.IntegrationMethod == "NULL"){
    file.remove(paste0(params.ProjectName,params.IntegrationMethod,"IntegratedUMAP.pdf"))
}

# ╔═════════════════════════════════╗
# ╠═ Plot Integrated TSNE (if run) ═╣
# ╚═════════════════════════════════╝
pdf(paste0(params.ProjectName,params.IntegrationMethod,"IntegratedTSNE.pdf"),width = 20, height = 15)
if (params.IntegrationMethod != "NULL"){
    message("Generating Integrated TSNE Plots")
    p1 <- DimPlot(object = SubsetSO, reduction = paste0("tsne.",params.IntegrationMethod,".Subset"), pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = "orig.ident", shuffle = T)
    p2 <- FeaturePlot(SubsetSO, reduction = paste0("tsne.",params.IntegrationMethod,".Subset"), pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "nFeature_RNA", order = T) + scale_color_viridis(limits =c(min(SubsetSO$nFeature_RNA),max(SubsetSO$nFeature_RNA)), direction = -1)
    p3 <- FeaturePlot(SubsetSO, reduction = paste0("tsne.",params.IntegrationMethod,".Subset"), pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "nCount_RNA", order = T) + scale_color_viridis(limits =c(min(SubsetSO$nCount_RNA),max(SubsetSO$nCount_RNA)), direction = -1)
    p4 <- FeaturePlot(SubsetSO, reduction = paste0("tsne.",params.IntegrationMethod,".Subset"), pt.size = (-0.00001837*length(SubsetSO$orig.ident))+1, features = "percent.mt", order = T) + scale_color_viridis(limits =c(min(SubsetSO$percent.mt),max(SubsetSO$percent.mt)), direction = -1)
    
    print(p1 + p2 + p3 + p4 + plot_layout(design = Page1layout))
    try(expr = {print(DimPlot(object = SubsetSO, reduction = paste0("tsne.",params.IntegrationMethod,".Subset"), pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = "CIscCATCH"))})
    for (i in params.Resolutions){
        print(DimPlot(object = SubsetSO, reduction = paste0("tsne.",params.IntegrationMethod,".Subset"), pt.size =(-0.00007653*length(SubsetSO$orig.ident))+4, label = T, group.by = paste0(params.IntegrationMethod,"SubsetRes.",i)))
    }
}
dev.off()

if(params.IntegrationMethod == "NULL"){
    file.remove(paste0(params.ProjectName,params.IntegrationMethod,"IntegratedTSNE.pdf"))
}

# ╔═══════════════════╗
# ╠═ Make Loupe File ═╣
# ╚═══════════════════╝
message("Generating Loupe File")
EULAmessage <- NULL
if(toupper(params.MakeLoupe) == "TRUE"){
    if (toupper(params.10xEULA) == "AGREE"){
        library(loupeR)
        loupeR::setup()
        create_loupe(count_mat = SubsetSO@assays$RNA$counts,
                     clusters = select_clusters(SubsetSO),
                     projections = select_projections(SubsetSO),
                     output_name = params.ProjectName
        )
    } else {
        EULAmessage <- "WARNING: Loupe File set to TRUE but you have not agreed to the 10x EULA -- Set params.10xEULA to Agree to create Loupe File."
        message(EULAmessage)
    }
}else {
    if (toupper(params.10xEULA) == "AGREE"){
        library(loupeR)
        loupeR::setup()
        create_loupe(count_mat = SubsetSO@assays$RNA$counts,
                     clusters = select_clusters(SubsetSO),
                     projections = select_projections(SubsetSO),
                     output_name = params.ProjectName
        )
        EULAmessage <- "NOTE: Loupe File set to FALSE but you have agreed to the 10x EULA -- Loupe File has been made."
        message(EULAmessage)
    }
}

sink(file = paste0(params.ProjectName,"SubsetExecution.log"), append = T)
cat("\n")
cat("╠═════════════════════════════════  Dimmensional Reduction  ═══════════════════════════════════╣\n")
cat(paste0("Final Seurat Object: \n"))
print(SubsetSO)
if (length(EULAmessage) != 0){
    cat(EULAmessage)
    cat("\n")
}
cat("\n")
cat("╠══════════════════════════════════════  Session Info  ════════════════════════════════════════╣\n")
sessionInfo()
sink()


