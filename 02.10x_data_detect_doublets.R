########01.We used DoubletFinder to detect doublets in each sample.
library(DoubletFinder)
####  set options 
args<-c()
args$percentage<- 0.05
args$dim<-20
doublets.percentage <- if(!is.null(args$percentage)) args$percentage else 0.05
dim.usage <- if(!is.null(args$dim)) args$dim else 20
args$out<-"./doublet"
#########################################
sample.all.qc=readRDS("sample.all.qc.rds")
table(sample.all.qc@meta.data$fetus)
batches<-unique(sample.all.qc@meta.data$fetus)

Find_doublet <- function(batch_i){
  batch_now<-batches[i]
  subset_seu<-subset(sample.all.qc,subset=sample==batch_now)
  subset_seu <- NormalizeData(subset_seu)
  subset_seu <- FindVariableFeatures(subset_seu, selection.method = "vst", nfeatures = 2000)
  subset_seu <- ScaleData(subset_seu)
  subset_seu <- RunPCA(subset_seu)
  ElbowPlot(subset_seu, ndims=20, reduction="pca") 
  subset_seu <- RunUMAP(subset_seu, dims = 1:20)
  sweep.res.list_kidney <- paramSweep_v3(subset_seu, PCs = 1:20, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn <- find.pK(sweep.stats_kidney)
  doublets.percentage=ncol(subset_seu)*8*1e-6
  nExp_poi <- round(doublets.percentage*nrow(subset_seu@meta.data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  subset_seu <- doubletFinder_v3(subset_seu, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(subset_seu@meta.data)[ncol(subset_seu@meta.data)] = "doublet_info"
  DimPlot(subset_seu, reduction = "umap", group.by = "doublet_info")
  table(subset_seu@meta.data$doublet_info.1)
  write.table(subset_seu@meta.data,paste0(args$out,"/batch_",batch_now,"_doublets_info.txt"),sep="\t",quote = F)
}
for (i in 1:length(batches)) {
  print(i)
  batch_r<-batches[i]
  print(batch_r)
  Find_doublet(batch_i = i)
}

########Remove the doublet cells in the FGC
sample.all.qc.list <- SplitObject(sample.all.qc, split.by = "fetus")
for (i in 1:length(sample.all.qc.list)) {
  sample.all.qc.list[[i]] <- NormalizeData(sample.all.qc.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  sample.all.qc.list[[i]] <- FindVariableFeatures(
    sample.all.qc.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
}
anchors1 <- FindIntegrationAnchors(object.list = sample.all.qc.list)
sample.all.qc.anchor <- IntegrateData(anchorset = anchors1)
DefaultAssay(sample.all.qc.anchor) <- "integrated"
sample.all.qc.anchor=ScaleData(sample.all.qc.anchor)
sample.all.qc.anchor <- RunPCA(sample.all.qc.anchor, npcs = 30, verbose = T)
DimPlot(sample.all.qc.anchor, reduction = "pca", group.by="fetus") 
ElbowPlot(sample.all.qc.anchor, ndims=25, reduction="pca") 

sample.all.qc.anchor <- FindNeighbors(sample.all.qc.anchor, reduction = "pca", dims = 1:20)
sample.all.qc.anchor <- FindClusters(sample.all.qc.anchor, resolution = 0.5)
sample.all.qc.anchor <- RunUMAP(sample.all.qc.anchor, reduction = "pca", dims = 1:20)

DimPlot(sample.all.qc.anchor, reduction = "umap", group.by = "fetus", pt.size = .1)
DimPlot(sample.all.qc.anchor, reduction = "umap", label = TRUE, pt.size = .1)
plot_grid(p1, p2)

DefaultAssay(sample.all.qc.anchor) <- "RNA"
sample.all.qc.anchor <- ScaleData(sample.all.qc.anchor)

diff.wilcox1 = FindAllMarkers(sample.all.qc.anchor)
markers1 = subset(diff.wilcox1,avg_logFC>0&p_val<0.05)
top10 = markers1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = subset(sample.all.qc.anchor,downsample=50), features = top10$gene,group.bar=T,size = 4)

f=c("POU5F1","DDX4","DAZL","MKI67","PCNA","ARX","DLK1","CYP17A1","STAR","INSL3","ACTA2","RGS5","AMH","SOX9","CDH5","EGFL7","CD14","PTPRC","NKG7")
fley=c("STAR","TCF21","NR2F1","MAFB","COL3A1","INHBA","PTCH1","INSL3","ARX","DLK1","PDGFRA","HSD3B2","CYP17A1","IGF1","IGF2","IGFBP5","CFD","IGFBP3","VIT")
fser=c("GATA3","GATA2","GATA4","ALX1","SOX9","AMH","ITGA6","WFDC2","BEX2","PRND","DMRT1","CITED1","DHH","HSD17B3")
fpmc=c("ACTA2","MYH11","TPM1","TPM2","TPM4","MYL9","PTCH1","PTCH2","GLI1","IGFBP6")
fec=c("VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN","PECAM1","CDH5")
fmc=c("CD4","CD14","CD163","S100A4","TYROBP","LYZ","RGS1")
ffgc=c("POU5F1","NANOG","TFAP2C","BLIMP1","NANOS3","DDX4","DAZL","PIWIL2","PIWIL4","SIX1")
for(i in c("f","fley","fser","fpmc","fec","fmc","ffgc")){
  DotPlot(sample.all.qc.anchor, features =c(HB.gene,get(i)) ) + RotatedAxis()
}
#By querying the genes found by FindAllMarkers Founction for each cluster and the expression of marker genes in various cell types in each clusters, it was determined that FGC belonged to groups 17 and 18
#############################find doublet in FGC
ksfgcmet <- subset(sample.all.qc.anchor@meta.data, seurat_clusters %in% c("17","18"))
ksfgc <- subset(sample.all.qc.anchor, cells=row.names(ksfgcmet))
ksfgc <- FindVariableFeatures(ksfgc, selection.method = "vst", nfeatures = 2000)
scale.genesfgc <- rownames(ksfgc)
ksfgc <- ScaleData(ksfgc, features = scale.genesfgc)
ksfgc <- RunPCA(ksfgc, features = VariableFeatures(ksfgc))
ElbowPlot(ksfgc, ndims=20, reduction="pca")
ksfgc <- FindNeighbors(ksfgc, dims = 1:15)
ksfgc <- FindClusters(ksfgc, resolution = 0.5)
addmargins(with(ksfgc@meta.data,table(seurat_clusters,fetus)))
ksfgc = RunUMAP(ksfgc, dims = 1:15)
DimPlot(ksfgc, reduction = "umap",label = F,group.by = "fetus")
DimPlot(ksfgc, reduction = "umap",label = TRUE)
DotPlot(ksfgc, features =c(HB.gene,f)) + RotatedAxis()
DotPlot(ksfgc, features =ffgc) + RotatedAxis()
FeaturePlot(ksfgc, features =c("POU5F1","DDX4","DLK1","AMH","EGFL7","CDH5")) + RotatedAxis()
delete.doublet2 <- WhichCells(ksfgc, expression = ARX > 1 | DLK1 > 1 |AMH > 1 |SOX9 >1 |STAR>1)  
saveRDS(delete.doublet2,file ="delete.doublet2.rds" )
saveRDS(sample.all.qc.anchor,file ="sample.all.qc.anchor.rds" )

########Remove the doublet cells found by DoubletFinder other than the FGC 
setwd("./doublet")
doublet.data <- rbindlist(lapply(list.files(), fread))
names(doublet.data)[1]<-"cell"
doublet=data.frame(cell=doublet.data$cell,pANN_0.25_0.3_348=doublet.data$pANN_0.25_0.3_348,doublet=doublet.data$doublet_info)
meta=sample.all.qc.anchor@meta.data
meta$cell=row.names(meta) 
meta=merge(meta,doublet,by="cell")
meta$doublet_info=NA
meta$doublet_info[meta$doublet == "Singlet"] <- "Singlet"
meta$doublet_info[meta$doublet == "Doublet" & meta$seurat_cluster %in% c("17","18")] <- "Singlet"
meta$doublet_info[meta$doublet == "Doublet" & meta$seurat_cluster != 17 & meta$seurat_cluster != 18] <- "Doublet"
table(meta$doublet_info)
delete.doublet1=row.names(meta[meta$doublet_info=="Doublet",])
saveRDS(delete.doublet1,file="delete.doublet1.rds")
delete.doublet=c(delete.doublet1,delete.doublet2)
saveRDS(delete.doublet,file="delete.doublet.rds")