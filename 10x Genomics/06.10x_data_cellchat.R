library(Seurat)
library(CellChat)
library(NMF)
library(reticulate)
## load seurat objects 
sample.finally=readRDS("sample.finally.rds")
sample.finally@meta.data$sex=factor(sample.finally@meta.data$sex,levels = c("Ctrl_M","KS"),ordered = T )

## Prepare required input data for CellChat analysis
scRNA.NC <- subset(sample.finally, subset = (sex%in%c("Ctrl_M")))
scRNA.KS <- subset(sample.finally, subset = (sex%in%c("KS")))
Idents(scRNA.NC)=scRNA.NC$celltype
Idents(scRNA.KS)=scRNA.KS$celltype
my_levels = c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","SPC","SC1","SC2","PMC","EC","MC","NK")
## Create a CellChat object
scRNA.NC <- createCellChat(scRNA.NC@assays$RNA@data, meta = scRNA.NC@meta.data, group.by = "celltype")
scRNA.NC@idents <- factor(scRNA.NC@idents, levels= my_levels)
levels(scRNA.NC@idents) 
scRNA.KS <- createCellChat(scRNA.KS@assays$RNA@data, meta = scRNA.KS@meta.data, group.by = "celltype")
scRNA.KS@idents <- factor(scRNA.KS@idents, levels= my_levels)
levels(scRNA.KS@idents) 

future::plan("multisession", workers = 10)  # do parallel
options(future.globals.maxSize = 8000 * 1024^2)

## running CellChat to infer cell-cell communication by KS data
cellchat <- scRNA.KS
## set the ligand-receptor interaction database 
cellchat@DB <- CellChatDB.human
## Identify over-expressed ligands or receptors.
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
## Infer cell-cell communication at a ligand-receptor pair level
cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = F)
## Filter the cell-cell communication
cellchat <- filterCommunication(cellchat, min.cells = 5)
## Infer cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
## Calculate aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
scRNA.KS <- cellchat

## running CellChat to infer cell-cell communication by male controls
cellchat <- scRNA.NC
## set the ligand-receptor interaction database 
cellchat@DB <- CellChatDB.human
## Identify over-expressed ligands or receptors.
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
## Infer cell-cell communication at a ligand-receptor pair level
cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = F)
## Filter the cell-cell communication
cellchat <- filterCommunication(cellchat, min.cells = 5)
## Infer cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
## Calculate aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
scRNA.NC <- cellchat

## merge CellChat object of each dataset together
scRNA.list <- list(Ctrl_M = scRNA.NC, KS = scRNA.KS)
cellchat <- mergeCellChat(scRNA.list, add.names = names(scRNA.list))
df.net <- subsetCommunication(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
saveRDS(scRNA.list,file="scRNA.list.rds")
saveRDS(cellchat,file="cellchat.rds")

## Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count", color.use = c("#00AFBB","#E7B800"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight", color.use = c("#00AFBB","#E7B800"))
gg1+gg2

## Set colors
color_use=c("#32CD32","#FFDA2B","#8EA1CC","#7570B9","#9370D8","#BA55D3","#6F3B9E","#00BFFF","#87CEFA","#4b74b2","#BC8F8F","#A9A9A9","#808080","#696969")
names(color_use)=c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","SPC","SC1","SC2","PMC","EC","MC","NK")

## Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = color_use)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.use = color_use)

## Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat,comparison = c(1,2),color.use = color_use)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(1,2),color.use = color_use)
p=gg1 + gg2

## Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(scRNA.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(scRNA.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(scRNA.list[[i]], title = names(scRNA.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

## Visually compare inferred cell-cell communication networks by circle plot
for (i in c("WNT","TGFb")) { 
  pathways.show <- c(i)
  weight.max <- getMaxWeight(scRNA.list, slot.name = c("netP"), attribute = pathways.show)
  par(mfrow = c(1,2), xpd=TRUE)
  for (j in 1:length(scRNA.list)) {
    netVisual_aggregate(scRNA.list[[j]],color.use = color_use, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(scRNA.list)[j]))
  }
}
## Visualize inferred communication networks of individual L-R pairs by circle plot
for (i in c("CDH2","SEMA7")) {
  pathways.show <- c(i)
  pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
  LR.show <- pairLR[1,] 
  weight.max <- getMaxWeight(scRNA.list, slot.name = c("net"), attribute = LR.show)
  par(mfrow = c(1,2), xpd=TRUE)
  plots=list()
  for (j in 1:length(scRNA.list)) {
    gg=netVisual_individual(scRNA.list[[j]], signaling = pathways.show , pairLR.use = LR.show, layout = "circle",color.use = color_use,edge.weight.max = weight.max[1], edge.width.max = 10,signaling.name = paste(pathways.show, names(scRNA.list)[j]))
    plots[[j]] <- gg
  }
}
## Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T, color.use = c("#00AFBB","#E7B800"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = T, color.use = c("#00AFBB","#E7B800"))
p <- gg1 + gg2
ggsave("ranknet.pdf", p, width = 12, height = 20, device = "pdf")
## Export results and p-values, and perform Benjamini-Hochberg multiple testing correction
ranknet=gg1[["data"]] 
ranknet$BH =p.adjust(ranknet$pvalues,method = "BH")
ranknet$p.adj.signif <- ifelse(ranknet$BH < 0.001, "***",
                            ifelse(ranknet$BH < 0.01, "**",
                                   ifelse(ranknet$BH < 0.05, "*", "ns")))
ranknet$p.signif <- ifelse(ranknet$pvalues < 0.001, "***",
                               ifelse(ranknet$pvalues < 0.01, "**",
                                      ifelse(ranknet$pvalues < 0.05, "*", "ns")))
ranknet$matches <- ranknet$p.adj.signif == ranknet$p.signif
write.table(ranknet, "ranknet_BH.xls", quote=F, sep="\t", col.names=T, row.names=F)
## Compare the overall information flow of each signaling pathway between FGC and Sertoli cell types
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = c(1,2,8,9,10), targets.use = c(1,2,8,9,10), stacked = T, do.stat = TRUE,color.use = c("#00AFBB","#E7B800"))
ggsave("ranknet_FGC_Ser.pdf",gg1, width = 6, height = 15, device = "pdf")
## Export results and p-values, and perform Benjamini-Hochberg multiple testing correction
ranknet=gg1[["data"]] 
ranknet$BH =p.adjust(ranknet$pvalues,method = "BH")
ranknet$p.adj.signif <- ifelse(ranknet$BH < 0.001, "***",
                            ifelse(ranknet$BH < 0.01, "**",
                                   ifelse(ranknet$BH < 0.05, "*", "ns")))
ranknet$p.signif <- ifelse(ranknet$pvalues < 0.001, "***",
                               ifelse(ranknet$pvalues < 0.01, "**",
                                      ifelse(ranknet$pvalues < 0.05, "*", "ns")))
ranknet$matches <- ranknet$p.adj.signif == ranknet$p.signif
write.table(ranknet, "ranknet_FGC_Ser_BH.xls", quote=F, sep="\t", col.names=T, row.names=F)

#####Compare the communication probabilities from certain cell groups to other cell groups
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CDH"))
pairLR.use  = pairLR.use[c(1),,drop=F]
a=netVisual_bubble(cellchat, comparison = c(1,2),sources.use = c(1,2,8,9,10), targets.use = c(1,2,8,9,10), 
                   pairLR.use = pairLR.use, remove.isolate = TRUE,max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, color.text = c("#00AFBB","#E7B800"))

pairLR.use <- extractEnrichedLR(cellchat, signaling = c("SEMA7"))
pairLR.use  = pairLR.use[c(1),,drop=F]
a=netVisual_bubble(cellchat, comparison = c(1,2),sources.use = c(1,2,8,9,10), targets.use = c(1,2,8,9,10), 
                   pairLR.use = pairLR.use, remove.isolate = TRUE,max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, color.text = c("#00AFBB","#E7B800"))

pairLR.use <- extractEnrichedLR(cellchat, signaling = c("JAM"))
pairLR.use  = pairLR.use[c(4),,drop=F]
a=netVisual_bubble(cellchat, comparison = c(1,2),sources.use = c(1,2,8,9,10), targets.use = c(1,2,8,9,10), 
                   pairLR.use = pairLR.use, remove.isolate = TRUE,max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, color.text = c("#00AFBB","#E7B800"))

pairLR.use <- extractEnrichedLR(cellchat, signaling = c("WNT"))
pairLR.use  = pairLR.use[c(81,86,5,45,13,53,82,87),,drop=F]
a=netVisual_bubble(cellchat, comparison = c(1,2), sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14) ,targets.use = c(1), 
                   pairLR.use = pairLR.use, remove.isolate = TRUE,max.dataset = 2, title.name = "Increased signaling", angle.x = 45, color.text = c("#00AFBB","#E7B800"))

pairLR.use <- extractEnrichedLR(cellchat, signaling = c("TGFb"))
pairLR.use  = pairLR.use[c(1,2,3,4,5,6,7,8,9),,drop=F]
a=netVisual_bubble(cellchat, comparison = c(1,2), sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14) ,targets.use = c(1), 
                   pairLR.use = pairLR.use,remove.isolate = TRUE,max.dataset = 2, title.name = "Increased signaling", angle.x = 45, color.text = c("#00AFBB","#E7B800"))
## bubble plot showing the Comparison result 
A=a[["data"]]
bk=c(0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30)
col_use<-c(colorRampPalette(colors = c('#FFE4E4',"#FF9595","#FF6060","#FF4646","#ff0000","#E30000"))(length(bk)/2),
           colorRampPalette(colors = c("#BE0000","#8b0000","#650000","#460000","#280000"))(length(bk)/2))
unique(A$group.names)
my_level=unique(A$group.names)
my_level <- my_level[length(my_level):1]
A$group.names=factor(A$group.names,levels = my_level,ordered = T)
print(
  ggplot(A, aes(x = dataset, y =group.names , size = pval, fill = prob, color = prob)) +
    geom_point(shape = 21, stroke = 1) +  
    scale_size_continuous(range = c(10),labels = c("p < 0.01"))+
    scale_fill_gradientn(colors = col_use) +  
    scale_color_gradientn(colors = col_use) +  
    labs(title = "KS Increased signaling", x = "Sex", y = "group") +  
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),  
          axis.text.x = element_text(size = 11),  
          axis.text.y = element_text(size = 11))+  
    facet_wrap(~interaction_name_2, ncol = 8)
)

print(
  ggplot(A, aes(x = dataset, y =group.names , size = pval, fill = prob, color = prob)) +
    geom_point(shape = 21, stroke = 1) +  
    scale_size_continuous(range = c(8, 10),breaks = c(2, 3),labels = c("0.01 < p < 0.05","p < 0.01"))+
    scale_fill_gradientn(colors = col_use) + 
    scale_color_gradientn(colors = col_use) + 
    labs(title = "KS Decreased signaling", x = "Sex", y = "group") + 
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),  
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11))+  
    facet_wrap(~interaction_name_2, ncol = 9)
)



