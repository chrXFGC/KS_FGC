library(readxl)
library(data.table) 
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)
library(Seurat)
library(Matrix)
library(ggpubr)
library(ggrepel)
library(scattermore)

################Create Seruat object and add the sex metadata
sample9 <- Read10X(data.dir = "filtered_feature_bc_matrix")
sample.all <- CreateSeuratObject(counts = sample9, project = "KS_10X", min.cells = 3, min.features = 200) 
sample_10x=data.frame(sample_name=c("KS3","KS4","KS6","Ctrl_M1","Ctrl_M2","Ctrl_M3","Ctrl_M4","Ctrl_M9","KS10"),row_num=c("1","2","6","3","4","5","7","8","9"))
for(i in 1:nrow(sample_10x)){
  sample.all@meta.data[which(substr(rownames(sample.all@meta.data), nchar(rownames(sample.all@meta.data)), nchar(rownames(sample.all@meta.data))) == sample_10x$row_num[i]),'sample'] <- sample_10x$sample_name[i]}
addmargins(table(sample.all@meta.data$sample)) 
sample.all@meta.data$sex <- ifelse(grepl("KS", sample.all@meta.data$sample), "KS", ifelse(grepl("Ctrl_M", sample.all@meta.data$sample), "Ctrl_M", ""))
addmargins(table(sample.all@meta.data$sex))

Ribo.gene=rownames(sample.all@assays[["RNA"]])[grepl("^RP[SL]", rownames(sample.all@assays[["RNA"]]))]#104
HB.gene=rownames(sample.all@assays[["RNA"]])[grepl("^HB[^(P)]", rownames(sample.all@assays[["RNA"]]))]#12
sample.all=PercentageFeatureSet(sample.all, pattern = "^MT-", col.name = "percent.mito")
sample.all=PercentageFeatureSet(sample.all, "^RP[SL]", col.name = "percent.ribo") 
sample.all=PercentageFeatureSet(sample.all, "^HB[^(P)]", col.name = "percent.hb")
saveRDS(sample.all,file ="sample.all.rds" )

###############Preliminary quality control
sample.all.qc <- subset(sample.all, subset = nFeature_RNA > 500 & percent.mito < 10 & percent.hb < 0.5)
addmargins(table(sample.all.qc@meta.data$sample)) 
addmargins(table(sample.all.qc@meta.data$sex))
saveRDS(sample.all.qc,file = "sample.all.qc.rds")

##############Analyze after removing doublet cells（delete.doublet can be accessed by running 02.10x_data_detect_doublets.R）
celllist=setdiff(colnames(sample.all.qc),delete.doublet) 
## Seruat object after removing doublet cells
sample.nod=sample.all[,celllist]
######## Integration of single-cell sequencing datasets based on Seurat package  ‘anchor-based’ integration workflow.
## split the object by dataset
sample.nod.list <- SplitObject(sample.nod, split.by = "sample")
## perform standard preprocessing on each object
for (i in 1:length(sample.nod.list)) {
  sample.nod.list[[i]] <- NormalizeData(sample.nod.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  sample.nod.list[[i]] <- FindVariableFeatures(
    sample.nod.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
}
## find anchors
anchors2 <- FindIntegrationAnchors(object.list = sample.nod.list)
## integrate data
sample.nod.anchor <- IntegrateData(anchorset = anchors2)
DefaultAssay(sample.nod.anchor) <- "integrated"
sample.nod.anchor=ScaleData(sample.nod.anchor)
sample.nod.anchor <- RunPCA(sample.nod.anchor, verbose = T)
DimPlot(sample.nod.anchor, reduction = "pca", group.by="sample") 
ElbowPlot(sample.nod.anchor, ndims=25, reduction="pca") 
sample.nod.anchor <- FindNeighbors(sample.nod.anchor, reduction = "pca", dims = 1:20)
sample.nod.anchor <- FindClusters(sample.nod.anchor, resolution = 0.5)
sample.nod.anchor <- RunUMAP(sample.nod.anchor, reduction = "pca", dims = 1:20)
DimPlot(sample.nod.anchor, reduction = "umap", group.by = "sample", pt.size = .1)
DimPlot(sample.nod.anchor, reduction = "umap", label = TRUE, pt.size = .1)

DefaultAssay(sample.nod.anchor) <- "RNA"
sample.nod.anchor <- ScaleData(sample.nod.anchor)

DotPlot(sample.nod.anchor, features =c(HB.gene,f) ) + RotatedAxis()
addmargins(with(sample.nod.anchor@meta.data,table(seurat_clusters,sample)))

## Find marker genes by FindAllMarkers Founction for each cluster
diff.wilcox2 = FindAllMarkers(sample.nod.anchor)
markers2 = subset(diff.wilcox2,avg_logFC>0&p_val<0.05)
top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = subset(sample.nod.anchor,downsample=50), features = top10$gene,group.bar=T,size = 4)

## Visualize the expression of marker genes for cell types in the testis across various clusters.
## "f","fley","fser","fpmc","fec","fmc","ffgc" defined in 02.10x_data_detect_doublets.R
for(i in c("f","fley","fser","fpmc","fec","fmc","ffgc")){
  DotPlot(sample.nod.anchor, features =c(HB.gene,get(i)) ) + RotatedAxis()
}
## By querying the genes identified by the FindAllMarkers function for each cluster and examining the expression of marker genes across different cell types within each cluster, clusters 20 and 21 were defined as undefined cells.
delete.undefined1=rownames(sample.nod.anchor@meta.data[sample.nod.anchor@meta.data$seurat_clusters %in% c(21,20),]) 
celllist=setdiff(celllist,delete.undefined1)
saveRDS(sample.nod,file = "sample.nod.rds")
saveRDS(sample.nod.anchor,file = "sample.nod.anchor.rds")

###########Analyze again after removing doublet cells and some undefined cells
## Seruat object after removing doublet cells and some undefined cells
sample.tmp=sample.all.qc[,celllist] 
sample.tmp.nohb<- sample.tmp[-grep("^HB[^(P)]",rownames(sample.tmp)),]
######## Integration of single-cell sequencing datasets based on Seurat package  ‘anchor-based’ integration workflow.
## split the object by dataset
sample.tmp.nohb.list <- SplitObject(sample.tmp.nohb, split.by = "sample")
## perform standard preprocessing on each object
for (i in 1:length(sample.tmp.nohb.list)) {
  sample.tmp.nohb.list[[i]] <- NormalizeData(sample.tmp.nohb.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  sample.tmp.nohb.list[[i]] <- FindVariableFeatures(
    sample.tmp.nohb.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
}
## find anchors
anchors.nohb <- FindIntegrationAnchors(object.list = sample.tmp.nohb.list)
## integrate data
sample.tmp.nohb_anchor <- IntegrateData(anchorset = anchors.nohb)
DefaultAssay(sample.tmp.nohb_anchor) <- "integrated"
sample.tmp.nohb_anchor=ScaleData(sample.tmp.nohb_anchor)
sample.tmp.nohb_anchor <- RunPCA(sample.tmp.nohb_anchor, verbose = T)
DimPlot(sample.tmp.nohb_anchor, reduction = "pca", group.by="sample") 
ElbowPlot(sample.tmp.nohb_anchor, ndims=40, reduction="pca") 
sample.tmp.nohb_anchor <- FindNeighbors(sample.tmp.nohb_anchor, reduction = "pca", dims = 1:20)
sample.tmp.nohb_anchor <- FindClusters(sample.tmp.nohb_anchor, resolution = 0.5)
sample.tmp.nohb_anchor <- RunUMAP(sample.tmp.nohb_anchor, reduction = "pca", dims = 1:20)
DimPlot(sample.tmp.nohb_anchor, reduction = "umap", label = TRUE, pt.size = .1)
DimPlot(sample.tmp.nohb_anchor, reduction = "umap", group.by = "sample", pt.size = .1)

DefaultAssay(sample.tmp.nohb_anchor) <- "RNA"
sample.tmp.nohb_anchor <- ScaleData(sample.tmp.nohb_anchor)

DotPlot(sample.tmp.nohb_anchor, features =c(f) ) + RotatedAxis()
addmargins(with(sample.tmp.nohb_anchor@meta.data,table(seurat_clusters,sample)))

## Find marker genes by FindAllMarkers Founction for each cluster
diff.wilcox3 = FindAllMarkers(sample.tmp.nohb_anchor)
markers3 = subset(diff.wilcox3,avg_logFC>0&p_val<0.05)
top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = subset(sample.tmp.nohb_anchor,downsample=50), features = top10$gene,group.bar=T,size = 4)

## Visualize the expression of marker genes for cell types in the testis across various clusters.
## "f","fley","fser","fpmc","fec","fmc","ffgc" defined in 02.10x_data_detect_doublets.R
for(i in c("f","fley","fser","fpmc","fec","fmc","ffgc")){
  DotPlot(sample.tmp.nohb_anchor, features =get(i) ) + RotatedAxis()
}
####################################subset leydig cells to reclustering
## By querying the genes identified by the FindAllMarkers function for each cluster and examining the expression of marker genes across different cell types within each cluster, clusters "0","1","2","3","5","6","7","8","14",and "12" were defined as leydig cells.
ksley2met <- subset(sample.tmp.nohb_anchor@meta.data, seurat_clusters %in% c("0","1","2","3","5","6","7","8","14","12"))
ksley2 <- subset(sample.tmp.nohb_anchor, cells=row.names(ksley2met))
ley2.list <- SplitObject(ksley2, split.by = "sample")
for (i in 1:length(ley2.list)) {
  ley2.list[[i]] <- NormalizeData(ley2.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  ley2.list[[i]] <- FindVariableFeatures(
    ley2.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
}
ley2.anchors <- FindIntegrationAnchors(object.list = ley2.list)
ksley2_anchor <- IntegrateData(anchorset = ley2.anchors)
DefaultAssay(ksley2_anchor) <- "integrated"
ksley2_anchor=ScaleData(ksley2_anchor)
ksley2_anchor <- RunPCA(ksley2_anchor, verbose = T)
DimPlot(ksley2_anchor, reduction = "pca", group.by="sample")
ElbowPlot(ksley2_anchor, ndims=25, reduction="pca")
ksley2_anchor <- FindNeighbors(ksley2_anchor, reduction = "pca", dims = 1:20)
ksley2_anchor <- FindClusters(ksley2_anchor, resolution = 0.5)
ksley2_anchor <- RunUMAP(ksley2_anchor, reduction = "pca", dims = 1:20)
DimPlot(ksley2_anchor, reduction = "umap", group.by = "sample", pt.size = .1)
DimPlot(ksley2_anchor, reduction = "umap", label = TRUE, pt.size = .1)

DefaultAssay(ksley2_anchor) <- "RNA"
ksley2_anchor <- ScaleData(ksley2_anchor)
CellstoHighlight <- WhichCells(ksley2_anchor, idents = "10")
DimPlot(ksley2_anchor, reduction = "umap", label = TRUE, repel = TRUE, cells.highlight = row.names(ksley2_anchor@meta.data[ksley2_anchor@meta.data$seurat_clusters=="10",]), cols.highlight = "red") + NoLegend()
DotPlot(ksley2_anchor, features =c("STAR","INSL3","CYP17A1","HSD3B2") ) + RotatedAxis()
DotPlot(ksley2_anchor, features =f) + RotatedAxis()
DotPlot(ksley2_anchor, features =fley) + RotatedAxis()
addmargins(with(ksley2_anchor@meta.data,table(seurat_clusters,sample)))
FeaturePlot(ksley2_anchor, features =c("DLK1","AMH","SOX9","NR2F1") ) + RotatedAxis()

diff.wilcox.ley2 = FindAllMarkers(ksley2_anchor)
ley2.markers = diff.wilcox.ley2 %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = ley2.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(object = subset(ksley2_anchor,downsample=50), features = top10$gene,group.bar=T,size = 2)

saveRDS(ksley2_anchor,file = "ksley2_anchor.rds")
#####################################subset sertoli cells to reclustering
## By querying the genes identified by the FindAllMarkers function for each cluster and examining the expression of marker genes across different cell types within each cluster, clusters "13","4",and "10" were defined as sertoli cells.
ksser2met <- subset(sample.tmp.nohb_anchor@meta.data, seurat_clusters %in% c("13","4","10"))
ksser2 <- subset(sample.tmp.nohb_anchor, cells=row.names(ksser2met))
ser2.list <- SplitObject(ksser2, split.by = "sample")
for (i in 1:length(ser2.list)) {
  ser2.list[[i]] <- NormalizeData(ser2.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  ser2.list[[i]] <- FindVariableFeatures(
    ser2.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
}
k.filter <- min(200, min(sapply(ser2.list, ncol)))
ser2.anchors <- FindIntegrationAnchors(object.list = ser2.list,k.filter = k.filter)
ksser2_anchor <- IntegrateData(anchorset = ser2.anchors)
DefaultAssay(ksser2_anchor) <- "integrated"
ksser2_anchor=ScaleData(ksser2_anchor)
ksser2_anchor <- RunPCA(ksser2_anchor, verbose = T)
DimPlot(ksser2_anchor, reduction = "pca", group.by="sample")
ElbowPlot(ksser2_anchor, ndims=25, reduction="pca")
ksser2_anchor <- FindNeighbors(ksser2_anchor, reduction = "pca", dims = 1:15)
ksser2_anchor <- FindClusters(ksser2_anchor, resolution = 0.5)
ksser2_anchor <- RunUMAP(ksser2_anchor, reduction = "pca", dims = 1:15)
DimPlot(ksser2_anchor, reduction = "umap", group.by = "sample", pt.size = .1)
DimPlot(ksser2_anchor, reduction = "umap", label = TRUE, pt.size = .1)

DefaultAssay(ksser2_anchor) <- "RNA"
ksser2_anchor <- ScaleData(ksser2_anchor)
addmargins(with(ksser2_anchor@meta.data,table(seurat_clusters,sample)))
DotPlot(ksser2_anchor, features =c("EGR3","JUN","NR4A1","S100A13","ENO1","BEX1","HOPX","DEFB119","CST9L") ) + RotatedAxis()
DotPlot(ksser2_anchor, features =f) + RotatedAxis()
DotPlot(ksser2_anchor, features =fser) + RotatedAxis()
FeaturePlot(ksser2_anchor, features =c("DLK1","AMH","XIST","MKI67")) + RotatedAxis()

diff.wilcox.ser2 = FindAllMarkers(ksser2_anchor)
ser2.markers = diff.wilcox.ser2 %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = ser2.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(object = subset(ksser2_anchor,downsample=50), features = top10$gene,group.bar=T,size = 2)
saveRDS(ksser2_anchor,file = "ksser2_anchor.rds")

## By querying the genes identified by the FindAllMarkers function for each cluster and examining the expression of marker genes across different cell types within each cluster, clusters "11" in leydig lineage and clusters "6" in sertoli lineage were defined as undefined cells.
delete.undefined2=row.names(ksser2_anchor@meta.data[ksser2_anchor@meta.data$seurat_clusters==6,])
delete.undefined3=row.names(ksley2_anchor@meta.data[ksley2_anchor@meta.data$seurat_clusters==11,])
delete.undefined=c(delete.undefined1,delete.undefined2,delete.undefined3)
#########################################subset FGC cells to reclustering
fgcmet <- subset(sample.tmp.nohb_anchor@meta.data, seurat_clusters %in% c("16"))
fgc <- subset(sample.tmp.nohb_anchor, cells=row.names(fgcmet))
table(fgc$sample)
fgc <- FindVariableFeatures(fgc, selection.method = "vst", nfeatures = 2000)
scale.genesfgc <- rownames(fgc)
fgc <- ScaleData(fgc, features = scale.genesfgc)
fgc <- RunPCA(fgc, features = VariableFeatures(fgc))
ElbowPlot(fgc, ndims=20, reduction="pca")
fgc <- FindNeighbors(fgc, dims = 1:13)
fgc <- FindClusters(fgc, resolution = 0.5)
addmargins(with(fgc@meta.data,table(seurat_clusters,sample)))
fgc = RunUMAP(fgc, dims = 1:13)
DimPlot(fgc, reduction = "umap",label = F,group.by = "sample")
DimPlot(fgc, reduction = "umap",label = TRUE)
DotPlot(fgc, features =f) + RotatedAxis()
DotPlot(fgc, features =ffgc) + RotatedAxis()
FeaturePlot(fgc, features =c("POU5F1","DDX4","DLK1","AMH","EGFL7","CDH5")) + RotatedAxis()

delete=c(delete.doublet,delete.undefined)
saveRDS(delete.undefined,file = "delete.undefined.rds")
saveRDS(delete,file = "delete.rds")
saveRDS(sample.tmp.nohb_anchor,file ="sample.tmp.nohb_anchor.rds" )
saveRDS(sample.tmp,file="sample.tmp.rds")

##########Analyze again after removing doublet cells and undefined cells
## Seruat object after removing doublet cells andundefined cells
celllist=setdiff(colnames(sample.all.qc),delete) 
sample.final=sample.all.qc[,celllist] 
saveRDS(sample.final,file = "sample.final.rds")
saveRDS(celllist,file = "celllist.rds")

sample.final=PercentageFeatureSet(sample.final, pattern = "^MT-", col.name = "percent.mito")
sample.final=PercentageFeatureSet(sample.final, "^RP[SL]", col.name = "percent.ribo") 
sample.final=PercentageFeatureSet(sample.final, "^HB[^(P)]", col.name = "percent.hb")
VlnPlot(sample.final, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo","percent.hb"), ncol = 5,pt.size = 0)

sample.final<- sample.final[-grep("^HB[^(P)]",rownames(sample.final)),]

######## Integration of single-cell sequencing datasets based on Seurat package  ‘anchor-based’ integration workflow.
## split the object by dataset
sample.final.list <- SplitObject(sample.final, split.by = "sample")
## perform standard preprocessing on each object
for (i in 1:length(sample.final.list)) {
  sample.final.list[[i]] <- NormalizeData(sample.final.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  sample.final.list[[i]] <- FindVariableFeatures(
    sample.final.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
}
## find anchors
anchors <- FindIntegrationAnchors(object.list = sample.final.list)
## integrate data
sample.final_anchor <- IntegrateData(anchorset = anchors)
DefaultAssay(sample.final_anchor) <- "integrated"
sample.final_anchor=ScaleData(sample.final_anchor)
sample.final_anchor <- RunPCA(sample.final_anchor, verbose = T)
DimPlot(sample.final_anchor, reduction = "pca", group.by="sample") 
ElbowPlot(sample.final_anchor, ndims=40, reduction="pca") 
sample.final_anchor <- FindNeighbors(sample.final_anchor, reduction = "pca", dims = 1:20)
sample.final_anchor <- FindClusters(sample.final_anchor, resolution = 0.5)
sample.final_anchor <- RunUMAP(sample.final_anchor, reduction = "pca", dims = 1:20)
DimPlot(sample.final_anchor, reduction = "umap", label = TRUE, pt.size = .1)
DimPlot(sample.final_anchor, reduction = "umap", group.by = "sample", pt.size = .1)
DefaultAssay(sample.final_anchor) <- "RNA"
sample.final_anchor <- ScaleData(sample.final_anchor)

DotPlot(sample.final_anchor, features =f ) + RotatedAxis()
addmargins(with(sample.final_anchor@meta.data,table(seurat_clusters,sample)))

## Find marker genes by FindAllMarkers Founction for each cluster
diff.wilcox3 = FindAllMarkers(sample.final_anchor)
markers3 = subset(diff.wilcox3,avg_logFC>0&p_val<0.05)
top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = subset(sample.final_anchor,downsample=50), features = top10$gene,group.bar=T,size = 4)

## Visualize the expression of marker genes for cell types in the testis across various clusters.
## "f","fley","fser","fpmc","fec","fmc","ffgc" defined in 02.10x_data_detect_doublets.R
for(i in c("f","fley","fser","fpmc","fec","fmc","ffgc")){
  DotPlot(sample.final_anchor, features =get(i) ) + RotatedAxis()
}

####################################subset leydig cells to reclustering
## By querying the genes identified by the FindAllMarkers function for each cluster and examining the expression of marker genes across different cell types within each cluster, clusters "0","1","2","3","5","6","7","11","12",and "14" were defined as leydig cells.
leymet <- subset(sample.final_anchor@meta.data, seurat_clusters %in% c("0","1","2","3","5","6","7","11","12","14"))
ley<- subset(sample.final_anchor, cells=row.names(leymet))
ley.list <- SplitObject(ley, split.by = "sample")
for (i in 1:length(ley.list)) {
  ley.list[[i]] <- NormalizeData(ley.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  ley.list[[i]] <- FindVariableFeatures(
    ley.list[[i]], selection.method = "vst",
    nfeatures = 4000, verbose = FALSE)
}
ley.anchors <- FindIntegrationAnchors(object.list = ley.list)
ley_anchor <- IntegrateData(anchorset = ley.anchors)
DefaultAssay(ley_anchor) <- "integrated"
ley_anchor=ScaleData(ley_anchor)
ley_anchor <- RunPCA(ley_anchor, verbose = T)
DimPlot(ley_anchor, reduction = "pca", group.by="sample")
ElbowPlot(ley_anchor, ndims=45, reduction="pca")
ley_anchor <- FindNeighbors(ley_anchor, reduction = "pca", dims = 1:16)
ley_anchor <- FindClusters(ley_anchor, resolution = 0.5)
ley_anchor <- RunUMAP(ley_anchor, reduction = "pca", dims = 1:16)
DimPlot(ley_anchor, reduction = "umap", group.by = "sample", pt.size = .1)
DimPlot(ley_anchor, reduction = "umap", label = TRUE, pt.size = .1)

DefaultAssay(ley_anchor) <- "RNA"
ley_anchor <- ScaleData(ley_anchor)
DotPlot(ley_anchor, features =c("STAR","INSL4","CYP17A1","HSD4B4") ) + RotatedAxis()
DotPlot(ley_anchor, features =f) + RotatedAxis()
DotPlot(ley_anchor, features =fley) + RotatedAxis()
addmargins(with(ley_anchor@meta.data,table(seurat_clusters,sample)))
FeaturePlot(ley_anchor, features =c("DLK1","AMH","ACTA2","NR2F1") ) + RotatedAxis()

## renaming the leydig cell clusters as: "LPC1","LPC2","LP1","LP2","DLC"
ley.celltype=data.frame(ClusterID=c(0,1,2,3,4,5,6,7,8,9,10),celltype=c("LP1","LP1","LP1","LP2","LP1","LP1","LP1","LPC2","LPC1","DLC","LP1"))
ley_anchor@meta.data$celltype = "NA"
for(i in 1:nrow(ley.celltype)){
  ley_anchor@meta.data[which(ley_anchor@meta.data$seurat_clusters == ley.celltype$ClusterID[i]),'celltype'] <- ley.celltype$celltype[i]}
ley_anchor@meta.data$celltype <- factor(ley_anchor@meta.data$celltype, levels = c("LPC1","LPC2","LP1","LP2","DLC"), ordered = TRUE)
DimPlot(ley_anchor, reduction = "umap", group.by = "celltype", pt.size = .1)
addmargins(with(ley_anchor@meta.data,table(celltype,sample)))

## Find marker genes by FindAllMarkers Founction for each leydig cell type
Idents(ley_anchor) <- ley_anchor@meta.data$celltype
levels(ley_anchor)
diff.wilcox.ley = FindAllMarkers(ley_anchor)
ley.markers = subset(diff.wilcox.ley,avg_logFC>0&p_val_adj<0.1)
ley.top10 = ley.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = subset(ley_anchor,downsample=50), features = ley.top10$gene,group.bar=T,size = 4)
write.table(ley.markers,file="Leydig.marker.xls",col.names=T,row.names=F,sep="\t",quote=F)

## Calculate the number and proportion of each cell type within the Leydig lineage
ley_anchor@meta.data$celltype=as.character(ley_anchor@meta.data$celltype)
ley.proportion <- ley_anchor@meta.data %>%
  group_by(sex,celltype) %>%
  tally() %>%
  transmute(sex,celltype,n,lineage_sum=sum(n),pct = 100*n / sum(n))
write.table(ley.proportion, file="leydig.celltype_rate.xls",col.names=T,row.names=F,sep="\t",quote=F)

## Visualize cell proportions of each cell type within the Leydig lineage
ley.proportion$celltype <- factor(ley.proportion$celltype, levels = c("LPC1","LPC2","LP1","LP2","DLC"), ordered = T)
ley.proportion$sex <- factor(ley.proportion$sex, levels = c("Ctrl_M","KS"), ordered = T)
ggbarplot(ley.proportion, x="sex", y="pct", fill = "celltype",color = "celltype",
          xlab="Group",ylab="Percent",legend="right")+
  ggtitle("sex")+
  rotate_x_text(45)+guides(fill=guide_legend(title="sex"))
saveRDS(ley_anchor,file = "ley_anchor.rds")
#####################################subset sertoli cells to reclustering
## By querying the genes identified by the FindAllMarkers function for each cluster and examining the expression of marker genes across different cell types within each cluster, clusters "4","10",and "13" were defined as sertoli cells
sermet <- subset(sample.final_anchor@meta.data, seurat_clusters %in% c("4","10","13"))
ser <- subset(sample.final_anchor, cells=row.names(sermet))
ser.list <- SplitObject(ser, split.by = "sample")
for (i in 1:length(ser.list)) {
  ser.list[[i]] <- NormalizeData(ser.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  ser.list[[i]] <- FindVariableFeatures(
    ser.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
}
k.filter <- min(200, min(sapply(ser.list, ncol)))
ser.anchors <- FindIntegrationAnchors(object.list = ser.list,k.filter = k.filter)
ser_anchor <- IntegrateData(anchorset = ser.anchors)
DefaultAssay(ser_anchor) <- "integrated"
ser_anchor=ScaleData(ser_anchor)
ser_anchor <- RunPCA(ser_anchor, verbose = T)
DimPlot(ser_anchor, reduction = "pca", group.by="sample")
ElbowPlot(ser_anchor, ndims=25, reduction="pca")
ser_anchor <- FindNeighbors(ser_anchor, reduction = "pca", dims = 1:14)
ser_anchor <- FindClusters(ser_anchor, resolution = 0.5)
ser_anchor <- RunUMAP(ser_anchor, reduction = "pca", dims = 1:14)
DimPlot(ser_anchor, reduction = "umap", group.by = "sample", pt.size = .1)
DimPlot(ser_anchor, reduction = "umap", label = TRUE, pt.size = .1)
DefaultAssay(ser_anchor) <- "RNA"
ser_anchor <- ScaleData(ser_anchor)
DotPlot(ser_anchor, features =c("EGR4","JUN","NR4A1","S100A14","ENO1","BEX1","HOPX","DEFB119","CST9L") ) + RotatedAxis()
DotPlot(ser_anchor, features =f) + RotatedAxis()
DotPlot(ser_anchor, features =fser) + RotatedAxis()
FeaturePlot(ser_anchor, features =c("ARX","DLK1","SOX9","AMH")) + RotatedAxis()
addmargins(with(ser_anchor@meta.data,table(seurat_clusters,sample)))

## renaming the sertoli cell clusters as: "SPC","SC1","SC2"
ser.celltype=data.frame(ClusterID=c(0,1,2,3,4,5,6,7),celltype=c("SC1","SC1","SC1","SC1","SC2","SC2","SPC","SPC"))
ser_anchor@meta.data$celltype = "NA"
for(i in 1:nrow(ser.celltype)){
  ser_anchor@meta.data[which(ser_anchor@meta.data$seurat_clusters == ser.celltype$ClusterID[i]),'celltype'] <- ser.celltype$celltype[i]}
ser_anchor@meta.data$celltype <- factor(ser_anchor@meta.data$celltype, levels = c("SPC","SC1","SC2"), ordered = TRUE)
DimPlot(ser_anchor, reduction = "umap", group.by = "celltype", pt.size = .1)
addmargins(with(ser_anchor@meta.data,table(celltype,sample)))

## Find marker genes by FindAllMarkers Founction for each sertoli cell type
Idents(ser_anchor) <- ser_anchor@meta.data$celltype
levels(ser_anchor)
diff.wilcox.ser = FindAllMarkers(ser_anchor)
ser.markers = subset(diff.wilcox.ser,avg_logFC>0&p_val_adj<0.1)
ser.top10 = ser.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = subset(ser_anchor,downsample=50), features = ser.top10$gene,group.bar=T,size = 4)
write.table(ser.markers,file="Sertoli.marker.xls",col.names=T,row.names=F,sep="\t",quote=F)

##########subset sertoli cells with Ctrl_M to reclustering
ser.nc.anchor=subset(ser_anchor,cell=row.names(ser_anchor@meta.data[ser_anchor@meta.data$sex=="Ctrl_M",]))
addmargins(table(ser.nc.anchor@meta.data$sample))
ser.nc.anchor <- SplitObject(ser.nc.anchor, split.by = "sample")
for (i in 1:length(ser.nc.anchor)) {
  ser.nc.anchor[[i]] <- NormalizeData(ser.nc.anchor[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  ser.nc.anchor[[i]] <- FindVariableFeatures(
    ser.nc.anchor[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
}
k.filter <- min(200, min(sapply(ser.nc.anchor, ncol)))
ser.nc.anchors <- FindIntegrationAnchors(object.list = ser.nc.anchor,k.filter = k.filter)
ser.nc.anchor <- IntegrateData(anchorset = ser.nc.anchors)
DefaultAssay(ser.nc.anchor) <- "integrated"
ser.nc.anchor=ScaleData(ser.nc.anchor)
DefaultAssay(ser.nc.anchor) <- "RNA"
ser.nc.anchor <- ScaleData(ser.nc.anchor)
## Find marker genes by FindAllMarkers Founction for each sertoli cell type with Ctrl_M
Idents(ser.nc.anchor) <- ser.nc.anchor@meta.data$celltype
levels(ser.nc.anchor)
diff.wilcox.ser.nc = FindAllMarkers(ser.nc.anchor)
qc.ser.nc.markers = subset(diff.wilcox.ser.nc,avg_logFC>0&p_val_adj<0.1)
write.table(qc.ser.nc.markers,file="Sertoli.marker.in.Ctrl_M.xls",col.names=T,row.names=F,sep="\t",quote=F)

saveRDS(ser_anchor,file = "ser_anchor.rds")

## Calculate the number and proportion of each cell type within the sertoli lineage
ser_anchor@meta.data$celltype=as.character(ser_anchor@meta.data$celltype)
ser.proportion <- ser_anchor@meta.data %>%
  group_by(sex,celltype) %>%
  tally() %>%
  transmute(sex,celltype,n,lineage_sum=sum(n),pct = 100*n / sum(n))
write.table(ser.proportion, file="sertoli.celltype_rate.xls",col.names=T,row.names=F,sep="\t",quote=F)
unique(ser.proportion$celltype)

## Visualize cell proportions of each cell type within the sertoli lineage
ser.proportion$celltype <- factor(ser.proportion$celltype, levels = c("SPC","SC1","SC2"), ordered = T)
ser.proportion$sex <- factor(ser.proportion$sex, levels = c("Ctrl_M","KS"), ordered = T)
ggbarplot(ser.proportion, x="sex", y="pct", fill = "celltype",color = "celltype",
          xlab="Group",ylab="Percent",legend="right")+
  ggtitle("sex")+
  rotate_x_text(45)+guides(fill=guide_legend(title="sex"))
#####################################subset FGC cells to reclustering
fgcmet <- subset(sample.final_anchor@meta.data, seurat_clusters %in% c("17"))
fgc <- subset(sample.final_anchor, cells=row.names(fgcmet))#453
fgc <- FindVariableFeatures(fgc, selection.method = "vst", nfeatures = 2000)
scale.genesfgc <- rownames(fgc)
fgc <- ScaleData(fgc, features = scale.genesfgc)
fgc <- RunPCA(fgc, features = VariableFeatures(fgc))
ElbowPlot(fgc, ndims=20, reduction="pca")
fgc <- FindNeighbors(fgc, dims = 1:13)
fgc <- FindClusters(fgc, resolution = 0.5)
addmargins(with(fgc@meta.data,table(seurat_clusters,sample)))
fgc = RunUMAP(fgc, dims = 1:13)
DimPlot(fgc, reduction = "umap",label = F,group.by = "sample")
DimPlot(fgc, reduction = "umap",label = TRUE)

DotPlot(fgc, features =f) + RotatedAxis()
DotPlot(fgc, features =ffgc) + RotatedAxis()
FeaturePlot(fgc, features =c("TFAP2C","POU5F1","DDX4","DAZL"),ncol = 2) + RotatedAxis()
FeaturePlot(fgc, features =c("ARX","DLK1","AMH","SOX9"),ncol = 2) + RotatedAxis()
FeaturePlot(fgc, features =c("STAR","CYP17A1"),ncol = 2) + RotatedAxis()
double <- WhichCells(fgc, expression = ARX > 1 | DLK1 > 1 |AMH > 1 |SOX9 >1 |STAR>1 )

fgc <- subset(fgc, cells=setdiff(row.names(fgcmet),double) )
addmargins(with(ksfgc@meta.data,table(seurat_clusters,double)))

## renaming the FGC cell clusters as: "Early_FGC","Late_FGC"
fgc.celltype=data.frame(ClusterID=c(0,1,2,3,4,5),celltype=c("Late_FGC","Early_FGC","Early_FGC","Early_FGC","Early_FGC","Late_FGC"))
fgc@meta.data$celltype = "NA"
for(i in 1:nrow(fgc.celltype)){
  fgc@meta.data[which(fgc@meta.data$seurat_clusters == fgc.celltype$ClusterID[i]),'celltype'] <- fgc.celltype$celltype[i]}
fgc@meta.data$celltype <- factor(fgc@meta.data$celltype, levels = c("Early_FGC","Late_FGC"), ordered = TRUE)
DimPlot(fgc, reduction = "umap", group.by = "celltype", pt.size = .1)
addmargins(with(fgc@meta.data,table(celltype,sample)))

## Find marker genes by FindAllMarkers Founction for each FGC cell type 
Idents(fgc) <- fgc@meta.data$celltype
levels(fgc)
diff.wilcox.fgc = FindAllMarkers(fgc)
fgc.markers = subset(diff.wilcox.fgc,avg_logFC>0&p_val_adj<0.1)
top10 = fgc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = subset(fgc,downsample=50), features = top10$gene,group.bar=T,size = 2)
write.table(fgc.markers,file="fgc.marker.xls",col.names=T,row.names=F,sep="\t",quote=F)

## Calculate the number and proportion of each cell type within the FGC lineage
fgc.proportion <- fgc@meta.data %>%
  group_by(sex,celltype) %>%
  tally() %>%
  transmute(sex,celltype,n,lineage_sum=sum(n),pct = 100*n / sum(n))
write.table(fgc.proportion, file="fgc.celltype_rate.xls",col.names=T,row.names=F,sep="\t",quote=F)

## Visualize cell proportions of each cell type within the FGC lineage
fgc.proportion$celltype <- factor(fgc.proportion$celltype, levels = c("Early_FGC","Late_FGC"), ordered = T)
fgc.proportion$sex <- factor(fgc.proportion$sex, levels = c("Ctrl_M","KS"), ordered = T)
ggbarplot(fgc.proportion, x="sex", y="pct", fill = "celltype",color = "celltype",
          xlab="Group",ylab="Percent",legend="right")+
  ggtitle("sex")+
  rotate_x_text(45)+guides(fill=guide_legend(title="sex"))
##############################################finally cell type
## renaming all cell clusters as: "Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"
sample.final_anchors=subset(sample.final_anchor, cells=setdiff(row.names(sample.final_anchor@meta.data),double))
metadata=sample.final_anchors@meta.data
metadata$cell=row.names(metadata)
ser.mate=data.frame(cell=row.names(ser_anchor@meta.data),celltype=ser_anchor@meta.data$celltype)
ley.mate=data.frame(cell=row.names(ley_anchor@meta.data),celltype=ley_anchor@meta.data$celltype)
fgc.mate=data.frame(cell=row.names(fgc@meta.data),celltype=fgc@meta.data$celltype)
slf.mate=rbind(sl.mate,fgc.mate)
metadata$celltype=NA
for(i in 1:nrow(slf.mate)){
  metadata[which(metadata$cell == slf.mate$cell[i]),'celltype'] <- as.character(slf.mate$celltype[i])}
metadata$celltype=ifelse(metadata$seurat_clusters==15,"PMC",
                            ifelse(metadata$seurat_clusters==9,"EC",
                                   ifelse(metadata$seurat_clusters %in% c(16,8),"MC",
                                          ifelse(metadata$seurat_clusters==18,"NK",metadata$celltype))))
metadata$celltype <- factor(metadata$celltype, levels = c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
sample.final_anchors@meta.data=metadata
DimPlot(sample.final_anchors, reduction = "umap", group.by = "celltype", pt.size = .1,label = T)
addmargins(with(sample.final_anchors@meta.data,table(celltype,sample)))

## Calculate the number and proportion of each cell type 
proportion <- sample.final_anchors@meta.data %>%
  group_by(sex,celltype) %>%
  tally() %>%
  transmute(sex,celltype,n,lineage_sum=sum(n),pct = 100*n / sum(n))
write.table(proportion, file="all.celltype_rate.xls",col.names=T,row.names=F,sep="\t",quote=F)
proportion$celltype <- factor(proportion$celltype, levels = c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = T)
proportion$sex <- factor(proportion$sex, levels = c("Ctrl_M","KS"), ordered = T)

Idents(sample.final_anchors) <- sample.final_anchors@meta.data$celltype
levels(sample.final_anchors)
saveRDS(sample.final_anchors,file = "sample.finally.rds")

#############View the data for the Ctrl_M and KS separately.
sample.final_anchors.ks=subset(sample.final_anchors, cells=row.names(sample.final_anchors@meta.data[sample.final_anchors@meta.data$sex=="KS",]))
sample.final_anchors.Ctrl_M=subset(sample.final_anchors, cells=row.names(sample.final_anchors@meta.data[sample.final_anchors@meta.data$sex=="Ctrl_M",]))
DimPlot(sample.final_anchors.ks, reduction = "umap", group.by = "celltype", pt.size = .1,label = T)
DimPlot(sample.final_anchors.Ctrl_M, reduction = "umap", group.by = "celltype", pt.size = .1,label = T)
save(sample.final_anchors.ks,file =  "sample.finally.ks.rds")
save(sample.final_anchors.Ctrl_M,file =  "sample.finally.Ctrl_M.rds")

##############plot
sample.finally=readRDS("sample.finally.rds")
## define colors
color_use=c("#32CD32","#FFDA2B","#8EA1CC","#7570B9","#9370D8","#BA55D3","#6F3B9E","#00BFFF","#87CEFA","#4b74b2","#BC8F8F","#A9A9A9","#808080","#696969")
names(color_use)=c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","SPC","SC1","SC2","PMC","EC","MC","NK")
sex_col<- c("#00AFBB","#E7B800")
names(sex_col)=c("Ctrl_M","KS")

sample.finally@meta.data$sex=factor(sample.finally@meta.data$sex,levels = c("Ctrl_M","KS"),ordered = T )
Idents(sample.finally)=factor(sample.finally$celltype,levels = c("NK","MC","EC","SC2","SC1","SPC","PMC","DLC","LP2","LP1","LPC2","LPC1","Late_FGC","Early_FGC"), ordered = TRUE)

## Dot plot of marker genes for different cell types
pdf(file = "marker.dotplot.pdf",height = 5,width =8 )
fmarker=c("POU5F1","DDX4","MKI67","PCNA","ARX","DLK1","CYP17A1","STAR","ACTA2","RGS5","AMH","SOX9","CDH5","EGFL7","CD14","PTPRC","NKG7")
DotPlot(sample.finally, features = fmarker, cols = c("grey", "darkorange")) + RotatedAxis()+
  ylab("Cell Type")+
  xlab("Gene")
dev.off()

## UMAP plot for different cell types facet by sex
Idents(sample.finally)=factor(sample.finally$celltype,levels = c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
umap_result_coord <- as.data.frame(sample.finally@reductions$umap@cell.embeddings)
umap_result_coord$Cluster <- sample.finally@meta.data[rownames(umap_result_coord),"celltype"]
plot_df_corrd<-cbind(sample.finally@meta.data,umap_result_coord)
pdf(file = "umap_10x_bysex.pdf",height = 5,width =10 )
ggscatter(plot_df_corrd,x="UMAP_1",y="UMAP_2",color =  "celltype",
          size = 0.1,facet.by = "sex")+ 
  scale_color_manual(breaks = c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), values=color.use) +
  theme(legend.position = "right", 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),strip.text = element_text(color = "black",size = 12))+ 
  scale_size_continuous(range = c(0.1, 3))+
  guides(color = guide_legend(override.aes = list(size = 6)))
dev.off()

## UMAP plot for different cell types
pdf(file = "umap_10xall.pdf",height = 5,width =5.5 )
ggscatter(plot_df_corrd, x = "UMAP_1", y = "UMAP_2", color = "celltype", palette = color_use,
          shape = 19, size = 0.1, alpha = 0.7) +
  guides(color = guide_legend(override.aes = list(size = 6))) +  
  scale_size(range = c(0.1, 0.1))+
  theme(legend.position = "right")+
  theme_void()+
  geom_segment(aes(x = -7, y = -16, xend = -2.5, yend = -16),
               arrow = arrow(type = "closed", length = unit(0.05, "inches")),
               color = "black", size = 0.5)+
  geom_segment(aes(x = -7, y = -16, xend = -7, yend = -10),
               arrow = arrow(type = "closed", length = unit(0.05, "inches")),
               color = "black", size = 0.5)+
  annotate("text", x = -5, y = -17.5, label = "UMAP_1", vjust = -0.5) +
  annotate("text", x = -7.5, y = -16.5, label = "UMAP_2", hjust = -0.5,srt = 90, adj = 1)
dev.off()

## Visualize cell proportions of each cell type 
pdf(file = "celltype_rate_sex_10x.pdf",height = 5,width =4 )
ggbarplot(proportion, x="sex", y="pct", fill = "celltype",color = "celltype",palette = color_use,
          xlab="Sex",ylab="Percent",legend="right")+
  rotate_x_text(45)
dev.off()

## Visualize the cell proportions of each cell type for each sample
proportion.fetus <- sample.finally@meta.data %>%
  group_by(sample,celltype) %>%
  tally() %>%
  transmute(sample,celltype,n,lineage_sum=sum(n),pct = 100*n / sum(n))
proportion.fetus$celltype<- factor(proportion.fetus$celltype, levels = c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = T)
pdf(file = "celltype_rate_fetus_10x.pdf",height = 5,width =7 )
ggbarplot(proportion.fetus, x="sample", y="pct", fill = "celltype",color = "celltype",palette = color_use,
          xlab="Fetus",ylab="Percent",legend="right")+
  rotate_x_text(45)
dev.off()
