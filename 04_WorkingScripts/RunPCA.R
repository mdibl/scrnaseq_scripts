#!/usr/local/bin/Rscript


# title: "SeuratV5_Normalize_QC.Rmd"
# author: "Ryan Seaman"
# date: "02/06/2024"

##################
# LOAD LIBRARIES #
##################

library(dplyr)
library(Matrix)
library(viridis)
library(tidyverse)
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(SeuratWrappers)
library(Seurat.utils)
library(SingleCellExperiment)
library(gprofiler2)
library(ggplot2)
library(ggsankey)
library(DoubletFinder)
library(stringr)
library(patchwork)
library(loupeR)
library(presto)


###################################
# READ IN PARAMS AND DIRECTORIES #
##################################

args <- commandArgs(trailingOnly = TRUE)

# Load RDS
rds <- LoadSeuratRds(args[1])

# Maximum PC to be used in dimensional reduciton/ find neighbors
params.pcMax <- args[2]

# Sample Name
params.project_name <- args[3]

assign(params.project_name, rds)

###########
# RUN PCA #
###########

all.genes <- rownames(get(params.project_name))
SO <- RunPCA(get(params.project_name), features = all.genes, npcs = 100)

###############
# FIND MAX PC #
###############

Elbow <- ElbowPlot(SO,  ndims = 100, reduction = "pca")
ElbowPoints <- Elbow$data

for (i in seq(0.1, 0.9, 0.1)){
  loess <- loess(stdev ~ dims,data=ElbowPoints, span = i)
  ElbowPoints[paste0("loessS", i)] <- loess$fitted
}

plot(stdev ~ dims,data=ElbowPoints,pch=19,cex=0.1)
colors <- c('red', 'orange', 'yellow', 'green', 'blue', 'purple', 'pink', 'turquoise', 'brown')
counter <- 1
j <- order(ElbowPoints$dims)
for (i in seq(0.1, 0.9, 0.1)){
  lines(ElbowPoints$dims[j],ElbowPoints[paste0("loessS",i)][j],col=colors[counter],lwd=1)
  counter <- counter + 1
}

for (i in seq(0.1, 0.9, 0.1)){
  first_deriv <- diff(ElbowPoints$stdev)/diff(ElbowPoints$dims)
  ElbowPoints[paste0("stddev_Der2_", i)][2:100] <- diff(first_deriv)/diff(ElbowPoints$dims)
}

idents <- seq(0.1, 0.9, 0.1)
SSE <- c()
deriv2_var <- c()
for (item in columns){
  if (length(grep("Der2", item)) > 0){
    variance <- var(Elbow_data[item], na.rm = T)
    deriv2_var <- append(deriv2_var, variance)
  }else if (length(grep("loess", item)) > 0 & !(length(grep("Der", item)) > 0)){
    err <- Elbow_data[item] - Elbow_data$stdev
    norms <- sum(((err - min(err))**2)/((max(err) - min(err))**2))
    SSE <- append(SSE, norms)
  }
}

scores <- SSE + (deriv2_var * 10000)
index <- which(scores == min(scores))
ident <- idents[index]

print(ident)

if (toupper(params.pcMax) == "NULL"){
  ### Identify the point where the Elbow Plot flattens out (Verify in Plot)
  print(Elbow)
  pcCount <- 1
  print(pcCount)
  while(std_dev_data[pcCount]-std_dev_data[pcCount+1]>0.01 | std_dev_data[pcCount+1]-std_dev_data[pcCount+2]>0.01){
    pcCount <- pcCount + 1
  }
  
  params.pcMax <- pcCount

}else{
  params.pcMax
}

pdf(paste0(params.project_name,"_Merged_ElbowPlot.pdf"), width = 20, height = 15)
ggplot(df, aes(x, y)) +
  geom_point() +
  geom_line(aes(x, std_dev_data)) +
  xlab("PC") +
  ylab("Std Dev") +
  geom_vline(xintercept = params.pcMax, linetype="dotted", 
             color = "red", size=1.5)
ggtitle("Polynomial Regression of Std Dev ~ PC")
dev.off()

############
# SAVE RDS #
############

SaveSeuratRds(get(params.project_name), file = paste0(params.project_name, "_SO.rds"))

sink(paste0(params.project_name,".validation.log"))

print(get(params.project_name))
cat("\n")
print(my_summary)
cat("\n")
cat(paste0("PC Max Selected at: ", params.pcMax))

sink()




