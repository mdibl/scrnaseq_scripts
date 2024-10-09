# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                       Title: Destiny.R                                     ═╣
# ╠═                                    Updated: Sept 30 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Run Destiny Trajectory Analysis                                                        ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝

# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(destiny)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(viridis)
library(patchwork)

# ╔═══════════════════╗
# ╠═ Create Log File ═╣
# ╚═══════════════════╝
sink(file = paste0(params.ProjectName,"DestinyExecution.log"), append = T)
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  Destiny Execution Log\n")
cat(paste0("╠  Destiny Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("╠═══════════════════════════════════  Parameter Settings  ═════════════════════════════════════╣\n")
cat(paste0("params.ProjectName: ", params.ProjectName,"\n"))
cat(paste0("params.seuratObject: \n"))
print(params.seuratObject)
cat(paste0("params.npcs: \n"))
cat(paste0("\t",params.npcs,"\n"))
cat(paste0("params.distancemetrics: \n"))
cat(paste0("\t",params.distancemetrics,"\n"))
sink()

# ╔══════════════════════╗
# ╠═ Create Directories ═╣
# ╚══════════════════════╝
message("Creating Marker Directories")
for (pc in params.npcs){
    dir.create((paste0("npcs_",pc)), showWarnings = F)
    for (distance in params.distancemetrics){
        dir.create((paste0("npcs_",pc,"/",distance)), showWarnings = F)
    }
}

# ╔═════════════════════╗
# ╠═ Convert SO to SCE ═╣
# ╚═════════════════════╝
DestinySCE <- as.SingleCellExperiment(params.seuratObject)
saveRDS(DestinySCE, paste0(params.ProjectName,"DestinySCE.rds"))

sink(file = paste0(params.ProjectName,"DestinyExecution.log"), append = T)
cat("╠═════════════════════════════════════  Convert to SCE  ═══════════════════════════════════════╣\n")
cat(paste0("Singe Cell Experiment Object: \n"))
print(DestinySCE)
sink()

# ╔═════════════════════════╗
# ╠═ Run Destiny Test Mode ═╣
# ╚═════════════════════════╝
for (pc in params.npcs){
    for(distance in params.distancemetrics){
        cat(paste0("Running Destiny with ",pc," PCs and ",distance," distance metric\n"))
        DestinyDM <- DiffusionMap(DestinySCE, verbose = T, n_pcs = pc, distance = distance)
        saveRDS(DestinyDM, paste0("npcs_",pc,"/",distance,"/DestinyDM_pc",pc,"_distance",distance,".rds"))
        DestinyDM_Frame <- data.frame(DC1 = eigenvectors(DestinyDM)[, 1],
                                      DC2 = eigenvectors(DestinyDM)[, 2],
                                      DC3 = eigenvectors(DestinyDM)[, 3],
                                      DC4 = eigenvectors(DestinyDM)[, 4],
                                      orig.ident = DestinySCE$orig.ident)
        p1 <- ggplot(DestinyDM_Frame, aes(x = DC1, y = DC2, color = orig.ident)) + geom_point() + theme_minimal()
        p2 <- ggplot(DestinyDM_Frame, aes(x = DC1, y = DC3, color = orig.ident)) + geom_point() + theme_minimal()
        p3 <- ggplot(DestinyDM_Frame, aes(x = DC1, y = DC4, color = orig.ident)) + geom_point() + theme_minimal()
        p4 <- ggplot(DestinyDM_Frame, aes(x = DC2, y = DC3, color = orig.ident)) + geom_point() + theme_minimal()
        p5 <- ggplot(DestinyDM_Frame, aes(x = DC2, y = DC4, color = orig.ident)) + geom_point() + theme_minimal()
        p6 <- ggplot(DestinyDM_Frame, aes(x = DC3, y = DC4, color = orig.ident)) + geom_point() + theme_minimal()

        pdf(paste0("npcs_",pc,"/",distance,"/DestinyDM_pc",pc,"_distance",distance,".pdf"), width = 20, height = 10)
        print(p1 + p2 + p3 + p4 + p5 + p6)
        dev.off()
    }
}

sink(file = paste0(params.ProjectName,"DestinyExecution.log"), append = T)
cat("╠══════════════════════════════════════  Session Info  ════════════════════════════════════════╣\n")
sessionInfo()
sink()