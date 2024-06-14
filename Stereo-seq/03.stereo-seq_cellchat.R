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
## load seurat objects obtained from running 02.stereo-seq_expression.R
c161=readRDS("Ctrl_M8_1.rds")
c162=readRDS("Ctrl_M8_2.rds")
c181=readRDS("Ctrl_M9_1.rds")
c182=readRDS("Ctrl_M9_2.rds")
C51=readRDS("KS12_1.rds")
c52=readRDS("KS12_2.rds")
c53=readRDS("KS12_3.rds")
c61=readRDS("KS13_1.rds")
c62=readRDS("KS13_2.rds")
## Extract the necessary information for CellChat from the male controls
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
## generate necessary Gene expression data
genes.common <- Reduce(intersect, list(rownames(data.in1), rownames(data.in2), rownames(data.in3), rownames(data.in4)))
data.input <- cbind(data.in1[genes.common, ], data.in2[genes.common, ],data.in3[genes.common, ],data.in4[genes.common, ])
## generate necessary assigned cell labels
meta <- rbind(meta1, meta2,meta3,meta4)
rownames(meta) <- colnames(data.input)
meta$labels <- factor(meta$labels, levels = levels(Idents(cor16_1)))
meta$slices <- factor(meta$slices, levels = c("cor161", "cor162", "cor181", "cor182"))
unique(meta$labels)
unique(meta$slices)
## Spatial locations of spots from stereo-seq
spatial.locs <- rbind(spatial.locs1, spatial.locs2, spatial.locs3, spatial.locs4)
rownames(spatial.locs) <- colnames(data.input)
spatial.locs=spatial.locs[,c(1,2)]
## Scale factors and spot diameters of stereo-seq
spatial.factors <- rbind(spatial.factors1, spatial.factors2, spatial.factors3, spatial.factors4)
rownames(spatial.factors) <- c("cor161", "cor162", "cor181", "cor182")
meta.nc=meta
## Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
## Set ligand-receptor interaction database
cellchat@DB <- CellChatDB.human
## Identify over-expressed ligands or receptors
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
## Infer cell-cell communication network at both ligand-receptor pair and signaling pathway levels
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)
## Filter the cell-cell communication
cellchat <- filterCommunication(cellchat, min.cells = 10)
## Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
## Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
spRNA.Ctrl_M <- cellchat

## Extract the necessary information for CellChat from the KS
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
## generate necessary Gene expression data
genes.common <- Reduce(intersect, list(rownames(data.in1), rownames(data.in2), rownames(data.in3), rownames(data.in4), rownames(data.in5)))
data.input <- cbind(data.in1[genes.common, ], data.in2[genes.common, ],data.in3[genes.common, ],data.in4[genes.common, ],data.in5[genes.common, ])
## generate necessary assigned cell labels
meta <- rbind(meta1, meta2,meta3,meta4,meta5)
rownames(meta) <- colnames(data.input)
meta$labels <- factor(meta$labels, levels = levels(Idents(cor51)))
meta$slices <- factor(meta$slices, levels = c("cor61", "cor62", "cor51", "cor52", "cor53"))
meta.ks=meta
unique(meta$labels)
unique(meta$slices)
## Spatial locations of spots from stereo-seq
spatial.locs <- rbind(spatial.locs1, spatial.locs2, spatial.locs3, spatial.locs4, spatial.locs5)
rownames(spatial.locs) <- colnames(data.input)
spatial.locs=spatial.locs[,c(1,2)]
## Scale factors and spot diameters of stereo-seq
spatial.factors <- rbind(spatial.factors1, spatial.factors2, spatial.factors3, spatial.factors4, spatial.factors5)
rownames(spatial.factors) <- c("cor61", "cor62", "cor51", "cor52", "cor53")
## Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
## Set ligand-receptor interaction database
cellchat@DB <- CellChatDB.human
## Identify over-expressed ligands or receptors
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
## Infer cell-cell communication network at both ligand-receptor pair and signaling pathway levels
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)
## Filter the cell-cell communication
cellchat <- filterCommunication(cellchat, min.cells = 10)
## Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
## Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
spRNA.KS <- cellchat

## merge CellChat object of each dataset together
scRNA.list <- list(Ctrl_M = spRNA.Ctrl_M, KS = spRNA.KS)
cellchat <- mergeCellChat(scRNA.list, add.names = names(scRNA.list))
df.net <- subsetCommunication(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")

saveRDS(scRNA.list,file="spRNA.list.rds")
saveRDS(cellchat,file="sp.cellchat.rds")

##### Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets
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

#####Compare the communication probabilities from certain cell groups to other cell groups
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CDH"))
pairLR.use  = pairLR.use[c(1),,drop=F]
a=netVisual_bubble(cellchat, comparison = c(1,2),sources.use = c(1,2,9,10,11), targets.use = c(1,2,9,10,11),
                   pairLR.use = pairLR.use, remove.isolate = TRUE,max.dataset = 2, title.name = "Decreased signaling", angle.x = 45, color.text = c("#00AFBB","#E7B800"))
A=a[["data"]]
A$dataset=factor(A$dataset,levels = c("Ctrl_M","KS"),ordered = T)
## set colors
bk=c(0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30)
col_use<-c(colorRampPalette(colors = c('#FFE4E4',"#FF9595","#FF6060","#FF4646","#ff0000","#E30000"))(length(bk)/2),
           colorRampPalette(colors = c("#BE0000","#8b0000","#650000","#460000","#280000"))(length(bk)/2))
## bubble plot showing the Comparison result 
pdf(file="sptial_CDH2_buble_plot.pdf",width = 5,height = 7.5)
print(
    ggplot(A, aes(x = dataset, y = group.names, size = pval, fill = prob, color = prob)) +
      geom_point(shape = 21, stroke = 1) +
      scale_size_continuous(
        range = c(8, 10),
        breaks = c(2, 3),
        labels = c("0.01 < p < 0.05", "p < 0.01")
      ) +
      scale_fill_gradientn(colors = col_use) +  
      scale_color_gradientn(colors = col_use) +
      labs(title = "KS Decreased signaling", x = "Sex", y = "group") +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11)
      ) +
      facet_wrap(~interaction_name_2, ncol = 9)
)
dev.off()
