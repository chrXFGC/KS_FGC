library(readxl)
library(data.table) 
library(pheatmap)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(Matrix)
library(ggpubr)
library(ggrepel)
library(ggalluvial)
library(ggsci)
library(tidyverse)
library(CellChat)
library(NMF)
library(reticulate)
c161=readRDS("Ctrl_M8_1.rds")
c162=readRDS("Ctrl_M8_2.rds")
c181=readRDS("Ctrl_M9_1.rds")
c182=readRDS("Ctrl_M9_2.rds")
C51=readRDS("KS12_1.rds")
c52=readRDS("KS12_2.rds")
c53=readRDS("KS12_3.rds")
c61=readRDS("KS13_1.rds")
c62=readRDS("KS13_2.rds")
for (i in 1:4) {
  spot.size = 0.22
  conversion.factor = spot.size/0.44
  if (i <= 2) {
    cor <- get(paste("c16", i, sep = ""))
    col_idx <- 160 + i
    data.in=Seurat::GetAssayData(cor, slot = "data", assay = "SCT")
    colnames(data.in) <- paste0("cor",col_idx,"_",colnames(data.in))
    assign(paste0("data.in", i), data.in)
    meta = data.frame(labels = Idents(cor), slices = paste0("cor",col_idx))
    assign(paste0("meta", i), meta)
    spatial.locs = Seurat::GetTissueCoordinates(cor, scale = NULL, cols = c("imagerow", "imagecol")) 
    assign(paste0("spatial.locs", i), spatial.locs)
    spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
    assign(paste0("spatial.factors", i), spatial.factors)
  } else {
    cor <- get(paste("c18", i - 2, sep = ""))
    col_idx <- 180 + i - 2
    data.in=Seurat::GetAssayData(cor, slot = "data", assay = "SCT")
    colnames(data.in) <- paste0("cor",col_idx,"_",colnames(data.in))
    assign(paste0("data.in", i), data.in)
    meta = data.frame(labels = Idents(cor), slices = paste0("cor",col_idx))
    assign(paste0("meta", i), meta)
    spatial.locs = Seurat::GetTissueCoordinates(cor, scale = NULL, cols = c("imagerow", "imagecol")) 
    assign(paste0("spatial.locs", i), spatial.locs)
    spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
    assign(paste0("spatial.factors", i), spatial.factors)
  }
}
genes.common <- Reduce(intersect, list(rownames(data.in1), rownames(data.in2), rownames(data.in3), rownames(data.in4)))
data.input <- cbind(data.in1[genes.common, ], data.in2[genes.common, ],data.in3[genes.common, ],data.in4[genes.common, ])
meta <- rbind(meta1, meta2,meta3,meta4)
rownames(meta) <- colnames(data.input)
meta$labels <- factor(meta$labels, levels = levels(Idents(cor16_1)))
meta$slices <- factor(meta$slices, levels = c("cor161", "cor162", "cor181", "cor182"))
unique(meta$labels)
unique(meta$slices)
spatial.locs <- rbind(spatial.locs1, spatial.locs2, spatial.locs3, spatial.locs4)
rownames(spatial.locs) <- colnames(data.input)
spatial.locs=spatial.locs[,c(1,2)]
spatial.factors <- rbind(spatial.factors1, spatial.factors2, spatial.factors3, spatial.factors4)
rownames(spatial.factors) <- c("cor161", "cor162", "cor181", "cor182")
meta.nc=meta
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
spRNA.Ctrl_M <- cellchat

for (i in 1:5) {
  spot.size = 0.22
  conversion.factor = spot.size/0.44
  if (i <= 2) {
    cor <- get(paste("c6", i, sep = ""))
    col_idx <- 60 + i
    data.in=Seurat::GetAssayData(cor, slot = "data", assay = "SCT")
    colnames(data.in) <- paste0("cor",col_idx,"_",colnames(data.in))
    assign(paste0("data.in", i), data.in)
    meta = data.frame(labels = Idents(cor), slices = paste0("cor",col_idx))
    assign(paste0("meta", i), meta)
    spatial.locs = Seurat::GetTissueCoordinates(cor, scale = NULL, cols = c("imagerow", "imagecol")) 
    assign(paste0("spatial.locs", i), spatial.locs)
    spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
    assign(paste0("spatial.factors", i), spatial.factors)
  } else {
    cor <- get(paste("c5", i - 2, sep = ""))
    col_idx <- 50 + i - 2
    data.in=Seurat::GetAssayData(cor, slot = "data", assay = "SCT")
    colnames(data.in) <- paste0("cor",col_idx,"_",colnames(data.in))
    assign(paste0("data.in", i), data.in)
    meta = data.frame(labels = Idents(cor), slices = paste0("cor",col_idx))
    assign(paste0("meta", i), meta)
    spatial.locs = cor@images[["ST_bin40"]]@coordinates
    assign(paste0("spatial.locs", i), spatial.locs)
    spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
    assign(paste0("spatial.factors", i), spatial.factors)
  }
}
genes.common <- Reduce(intersect, list(rownames(data.in1), rownames(data.in2), rownames(data.in3), rownames(data.in4), rownames(data.in5)))
data.input <- cbind(data.in1[genes.common, ], data.in2[genes.common, ],data.in3[genes.common, ],data.in4[genes.common, ],data.in5[genes.common, ])
meta <- rbind(meta1, meta2,meta3,meta4,meta5)
rownames(meta) <- colnames(data.input)
meta$labels <- factor(meta$labels, levels = levels(Idents(cor51)))
meta$slices <- factor(meta$slices, levels = c("cor61", "cor62", "cor51", "cor52", "cor53"))
meta.ks=meta
unique(meta$labels)
unique(meta$slices)
spatial.locs <- rbind(spatial.locs1, spatial.locs2, spatial.locs3, spatial.locs4, spatial.locs5)
rownames(spatial.locs) <- colnames(data.input)
spatial.locs=spatial.locs[,c(1,2)]
spatial.factors <- rbind(spatial.factors1, spatial.factors2, spatial.factors3, spatial.factors4, spatial.factors5)
rownames(spatial.factors) <- c("cor61", "cor62", "cor51", "cor52", "cor53")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
spRNA.KS <- cellchat

scRNA.list <- list(Ctrl_M = spRNA.Ctrl_M, KS = spRNA.KS)
cellchat <- mergeCellChat(scRNA.list, add.names = names(scRNA.list))
df.net <- subsetCommunication(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")

saveRDS(scRNA.list,file="spRNA.list.rds")
saveRDS(cellchat,file="sp.cellchat.rds")

pathways.show <- "CDH"
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
LR.show <- pairLR[1,] # show one ligand-receptor pair
weight.max <- getMaxWeight(scRNA.list, slot.name = c("net"), attribute = "CDH2_CDH2")
netVisual_individual(scRNA.list[[1]], signaling = pathways.show , pairLR.use = "CDH2_CDH2", layout = "circle",color.use = color_use,edge.weight.max = weight.max[1], edge.width.max = 10,signaling.name = paste(pathways.show, names(scRNA.list)[1]))
netVisual_individual(scRNA.list[[2]], signaling = pathways.show , pairLR.use = "CDH2_CDH2", layout = "circle",color.use = color_use,edge.weight.max = weight.max[1], edge.width.max = 10,signaling.name = paste(pathways.show, names(scRNA.list)[2]))

pathways.show <- "SEMA7"
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
LR.show <- pairLR[1,] # show one ligand-receptor pair
weight.max <- getMaxWeight(scRNA.list, slot.name = c("net"), attribute = "SEMA7A_ITGB1_ITGA1")
netVisual_individual(scRNA.list[[1]], signaling = pathways.show , pairLR.use = "SEMA7A_ITGB1_ITGA1", layout = "circle",color.use = color_use,edge.weight.max = weight.max[1], edge.width.max = 10,signaling.name = paste(pathways.show, names(scRNA.list)[1]))
netVisual_individual(scRNA.list[[2]], signaling = pathways.show , pairLR.use = "SEMA7A_ITGB1_ITGA1", layout = "circle",color.use = color_use,edge.weight.max = weight.max[1], edge.width.max = 10,signaling.name = paste(pathways.show, names(scRNA.list)[2]))

pathways.show <- "JAM"
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
LR.show <- pairLR[4,] # show one ligand-receptor pair
weight.max <- getMaxWeight(scRNA.list, slot.name = c("net"), attribute = "JAM2_JAM3")
netVisual_individual(scRNA.list[[1]], signaling = pathways.show , pairLR.use = "JAM2_JAM3", layout = "circle",color.use = color_use,edge.weight.max = weight.max[1], edge.width.max = 10,signaling.name = paste(pathways.show, names(scRNA.list)[1]))
netVisual_individual(scRNA.list[[2]], signaling = pathways.show , pairLR.use = "JAM2_JAM3", layout = "circle",color.use = color_use,edge.weight.max = weight.max[1], edge.width.max = 10,signaling.name = paste(pathways.show, names(scRNA.list)[2]))
