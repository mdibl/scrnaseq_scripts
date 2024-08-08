#!/opt/conda/envs/dolphinnext/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                           Title: FindNeighborsClustersMarkers.R                            ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes in the Seurat Object generated from Integration.R (if-run, else RunPCA.R). Runs  ═╣
# ╠═     find neighbors, find clusters, find markers, and builds a Sankey Plot.                 ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝
# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(dplyr)
library(stringr)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggsankey)
library(patchwork)
library(presto)

# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
args <- commandArgs(trailingOnly = TRUE)

# RDS file from QC
params.SeuratObject <- args[1]

# Resolutions
params.Resolutions <- args[2]
params.Resolutions <- as.numeric(unlist(strsplit(params.Resolutions, ",")))
print(params.Resolutions)

# PC Max
params.pcMax <- as.integer(args[3])

# Integration Method: Options( CCA, RPCA, Harmony, FastMNN, NULL) where NULL is to not run
params.IntegrationMethod <- args[4]

# Project Name
params.ProjectName <- args[5]

# Scale Options ( SD or default SCT)
params.scaleMethod <- args[6]


# ╔══════════════════════╗
# ╠═ Create Directories ═╣
# ╚══════════════════════╝
dir.create("markers")
dir.create("markers/unintegrated/")
if (params.IntegrationMethod != "NULL"){
    dir.create("markers/integrated/")
}

# ╔════════════════════╗
# ╠═ Load Seurat .rds ═╣
# ╚════════════════════╝
MergedSO <- readRDS(params.SeuratObject)


# ╔════════════════════════════════════╗
# ╠═ Find Neighbors and Find Clusters ═╣
# ╚════════════════════════════════════╝
# Generate Resolution Names -- unintegrated
count <- 1
countMax <- length(params.Resolutions)
params.UnintegratedResolutionsNames <- c()
while (count <= countMax){
    resVal <- paste0("unintegratedRes.",params.Resolutions[count])
    params.UnintegratedResolutionsNames <- append(params.UnintegratedResolutionsNames, resVal)
    count <- count + 1
}
params.UnintegratedResolutionsNames

# Find Neighbors -- unintegrated
MergedSO <- FindNeighbors(MergedSO, dims = 1:params.pcMax, reduction = "pca")

# Find Clusters -- unintegrated
MergedSO <- FindClusters(MergedSO, resolution = params.Resolutions, cluster.name = params.UnintegratedResolutionsNames)

if (params.IntegrationMethod != "NULL"){
    # Generate Resolution Names -- integrated
    count <- 1
    countMax <- length(params.Resolutions)
    params.IntegratedResolutionsNames <- c()
    while (count <= countMax){
        resVal <- paste0(params.IntegrationMethod,"Res.",params.Resolutions[count])
        params.IntegratedResolutionsNames <- append(params.IntegratedResolutionsNames, resVal)
        count <- count + 1
    }
    params.IntegratedResolutionsNames
    # Find Neighbors -- integrated
    MergedSO <- FindNeighbors(MergedSO, dims = 1:params.pcMax, reduction = paste0("integrated.",params.IntegrationMethod))
    # Find Cluster -- integrated
    MergedSO <- FindClusters(MergedSO, resolution = params.Resolutions, cluster.name = params.IntegratedResolutionsNames)
}

# ╔════════════════════════╗
# ╠═ Generate Sankey Plot ═╣
# ╚════════════════════════╝
mergedSOmetaRes <- MergedSO@meta.data
mergedSOmetaRes <- mergedSOmetaRes[,grepl("^unintegratedRes",names(mergedSOmetaRes))]
names(mergedSOmetaRes) <- gsub("unintegratedRes","res",names(mergedSOmetaRes))

df <- mergedSOmetaRes%>%
    make_long(names(mergedSOmetaRes))

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
    labs(title = paste0("Sankey: ",params.ProjectName, " Unintegrated Resolutions")) + 
    labs(caption =NULL ) + 
    labs(fill = 'Nodes')

### Save Sankey Plot
pdf(paste0(params.ProjectName,"UnintegratedSankeyPlot.pdf"), width = 20, height = 15)
pl
dev.off()


if (params.IntegrationMethod != "NULL"){
    mergedSOmetaRes <- MergedSO@meta.data
    pattern <- paste0(params.IntegrationMethod,"Res")
    mergedSOmetaRes <- mergedSOmetaRes[,grepl(pattern,names(mergedSOmetaRes))]
    names(mergedSOmetaRes) <- gsub(pattern,"res",names(mergedSOmetaRes))
    
    df <- mergedSOmetaRes%>%
        make_long(names(mergedSOmetaRes))
    
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
        labs(title = paste0("Sankey: ",params.ProjectName," ", params.IntegrationMethod," Integrated Resolutions")) + 
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

MergedSO[["RNA"]] <- JoinLayers(MergedSO[["RNA"]])

# ╔════════════════════════════════════╗
# ╠═ Find Markers and Save TSV output ═╣
# ╚════════════════════════════════════╝
if (params.scaleMethod != "SD"){
    MergedSO <- PrepSCTFindMarkers(MergedSO, assay = "SCT", verbose = TRUE)
}

for(i in 1:length(params.Resolutions)){
    print(paste0("Finding Unintegrated Markers at ",params.Resolutions[i]," resolution..."))
    Idents(MergedSO) <- paste0("unintegratedRes.",params.Resolutions[i])
    assign(paste0("UnintegratedMarkersRes",params.Resolutions[i]),FindAllMarkers(MergedSO, only.pos = T, logfc.threshold = 0.3))
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
        Idents(MergedSO) <- paste0(params.IntegrationMethod,"Res.",params.Resolutions[i])
        assign(paste0("IntegratedMarkersRes",params.Resolutions[i]),FindAllMarkers(MergedSO, only.pos = T, logfc.threshold = 0.3))
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

# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
SaveSeuratRds(MergedSO, file = paste0("06_",params.ProjectName, "_ClusterSO.rds"))

# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
sink(paste0("06_",params.ProjectName,"_ClusterValidation.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  FindNeighborsClustersMarkers.R log\n")
cat(paste0("╠  Analysis Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
cat("Seurat Object Status:\n")
print(MergedSO)
sink()

sink(paste0("06_",params.ProjectName,"_ClusterVersions.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  FindNeighborsClustersMarkers.R Versions"\n)
cat(paste0("╠  Analysis Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
sessionInfo()
sink()
