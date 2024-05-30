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
library(smerc)
library(findPC)

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
  lines(ElbowPoints$dims[j],ElbowPoints[,paste0("loessS",i)][j],col=colors[counter],lwd=1)
  counter <- counter + 1
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
  plot(SSE, deriv2_var)
  index <- which(SSE == min(SSE))
  ident <- paste0('loessS', idents[index])
}else {
  index <- 1
  ident <- paste0('loessS', idents[index])
}

if (toupper(params.pcMax) == "NULL"){
  ### Identify the point where the Elbow Plot flattens out (Verify in Plot)
  print(Elbow)
  pcCount <- 1
  print(pcCount)
  while(ElbowPoints[[ident]][pcCount]-ElbowPoints[[ident]][pcCount+1]>0.05 | ElbowPoints[[ident]][pcCount+1]-ElbowPoints[[ident]][pcCount+2]>0.05){
    pcCount <- pcCount + 1
  }
  
  params.pcMax <- pcCount

}else{
  params.pcMax
}

x = ElbowPoints$dims
y = ElbowPoints$stdev
df <- data.frame(x,y)

pc_tbl <- findPC(sdev = ElbowPoints[[ident]], number = 100, method = "all", figure = T)
params.pcMax <- mean(x = c(pc_tbl[1,2], pc_tbl[1,3], pc_tbl[1,4]))

pdf(paste0(params.project_name,"_Merged_ElbowPlot.pdf"), width = 20, height = 15)
ggplot(df, aes(x, y)) +
  geom_point() +
  geom_line(aes(x = seq(1,100), y = ElbowPoints[[ident]]), color = "green") +
  xlab("PC") +
  ylab("Std Dev") +
  geom_vline(xintercept = params.pcMax, linetype="dotted", 
             color = "red", size=1.5) +
  ggtitle(paste0("Loess Regression of Std Dev ~ PC  :  PC Chosen = ", params.pcMax))
dev.off()

############
# SAVE RDS #
############

SaveSeuratRds(get(params.project_name), file = paste0(params.project_name, "_SO.rds"))

loess <- loess(stdev ~ dims,data=ElbowPoints, span = idents[index])
my_summary <- summary(loess)

sink(paste0(params.project_name,".validation.log"))

print(get(params.project_name))
cat("\n")
print(my_summary)
cat("\n")
cat(paste0("First Derivative : ", pc_tbl[1,2], "\n"))
cat(paste0("Second Derivative : ", pc_tbl[1,3], "\n"))
cat(paste0("Preceding Residual : ", pc_tbl[1,4], "\n"))
cat("\n")
cat(paste0("PC Max Selected at: ", params.pcMax))

sink()




