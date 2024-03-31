library(Seurat)
library(CellChat)
library(NMF)
library(reticulate)

sample.finally=readRDS("sample.finally.rds")
sample.finally@meta.data$sex=factor(sample.finally@meta.data$sex,levels = c("Ctrl_M","KS"),ordered = T )

scRNA.NC <- subset(sample.finally, subset = (sex%in%c("Ctrl_M")))
scRNA.KS <- subset(sample.finally, subset = (sex%in%c("KS")))
Idents(scRNA.NC)=scRNA.NC$celltype
Idents(scRNA.KS)=scRNA.KS$celltype
my_levels = c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","SPC","SC1","SC2","PMC","EC","MC","NK")
scRNA.NC <- createCellChat(scRNA.NC@assays$RNA@data, meta = scRNA.NC@meta.data, group.by = "celltype")
scRNA.NC@idents <- factor(scRNA.NC@idents, levels= my_levels)
levels(scRNA.NC@idents) 
scRNA.KS <- createCellChat(scRNA.KS@assays$RNA@data, meta = scRNA.KS@meta.data, group.by = "celltype")
scRNA.KS@idents <- factor(scRNA.KS@idents, levels= my_levels)
levels(scRNA.KS@idents) 

future::plan("multisession", workers = 10) 
options(future.globals.maxSize = 8000 * 1024^2)

cellchat <- scRNA.KS
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = F)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
scRNA.KS <- cellchat

cellchat <- scRNA.NC
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = F)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
scRNA.NC <- cellchat

scRNA.list <- list(Ctrl_M = scRNA.NC, KS = scRNA.KS)
cellchat <- mergeCellChat(scRNA.list, add.names = names(scRNA.list))
df.net <- subsetCommunication(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
saveRDS(scRNA.list,file="scRNA.list.rds")
saveRDS(cellchat,file="cellchat.rds")

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count", color.use = c("#00AFBB","#E7B800"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight", color.use = c("#00AFBB","#E7B800"))
gg1+gg2

color_use=c("#32CD32","#FFDA2B","#8EA1CC","#7570B9","#9370D8","#BA55D3","#6F3B9E","#00BFFF","#87CEFA","#4b74b2","#BC8F8F","#A9A9A9","#808080","#696969")
names(color_use)=c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","SPC","SC1","SC2","PMC","EC","MC","NK")

par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = color_use)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.use = color_use)

gg1 <- netVisual_heatmap(cellchat,comparison = c(1,2),color.use = color_use)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(1,2),color.use = color_use)
p=gg1 + gg2

num.link <- sapply(scRNA.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(scRNA.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(scRNA.list[[i]], title = names(scRNA.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

for (i in c("WNT","TGFb")) { 
  pathways.show <- c(i)
  weight.max <- getMaxWeight(scRNA.list, slot.name = c("netP"), attribute = pathways.show)
  par(mfrow = c(1,2), xpd=TRUE)
  for (j in 1:length(scRNA.list)) {
    netVisual_aggregate(scRNA.list[[j]],color.use = color_use, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(scRNA.list)[j]))
  }
}
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

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T, color.use = c("#00AFBB","#E7B800"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = T, color.use = c("#00AFBB","#E7B800"))
p <- gg1 + gg2
ggsave("ranknet.pdf", p, width = 12, height = 20, device = "pdf")
write.table(gg1[["data"]], "ranknet.xls", quote=F, sep="\t", col.names=T, row.names=F)


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

A=a[["data"]]
bk=c(0,0.05,0.1,0.15,0.2,0.25,0.3)
col_use<-c(colorRampPalette(colors = c('darkblue',"lightblue","#E6E6FA","#FFF0F5"))(length(bk)/2),
           colorRampPalette(colors = c("#FFE4E1","indianred1","darkred"))(length(bk)/2))
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



