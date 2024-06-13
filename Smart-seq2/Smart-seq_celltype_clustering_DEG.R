# Clear environment
rm(list=ls())
.libPaths("/software/R_lib")

Packages <- c("dplyr", "ggplot2", "tibble","data.table","ggpubr","Seurat","ggsci","cowplot")
lapply(Packages, require, character.only = TRUE)

############# Create Seruat object ######
TPM<-read.table("/Data/FGC_TPM.xls")
metadata<-fread("/Data/FGC_Metadata.xls")%>%column_to_rownames(var = "sample_name")
nor_tpm<- apply(TPM, 2, function(x) log2(x+1))


KS_object<-CreateSeuratObject(nor_tpm, meta.data = metadata, min.cells = 10,project = "KS_scRNAseq")

sample.all.qc.list <- SplitObject(object = KS_object, split.by = "source")
for (i in 1:length(sample.all.qc.list)) {
  sample.all.qc.list[[i]] <- FindVariableFeatures(
    sample.all.qc.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
}

# find anchors
anchors1 <- FindIntegrationAnchors(object.list = sample.all.qc.list)

# integrate data
sample.all.qc.anchor <- IntegrateData(anchorset = anchors1)

DefaultAssay(sample.all.qc.anchor) <- "integrated"

sample.all.qc.anchor=ScaleData(sample.all.qc.anchor)

sample.all.qc.anchor <- RunPCA(sample.all.qc.anchor, npcs = 30, verbose = T)
sample.all.qc.anchor <- RunUMAP(sample.all.qc.anchor, reduction = "pca", dims = 1:20)
# Clustering

sample.all.qc.anchor <- FindNeighbors(sample.all.qc.anchor, reduction = "pca", dims = 1:20)
sample.all.qc.anchor <- FindClusters(sample.all.qc.anchor, resolution = 0.25)


p1<-DimPlot(sample.all.qc.anchor, reduction = "umap", group.by = "sex_type", pt.size = .1)
p2<-DimPlot(sample.all.qc.anchor, reduction = "umap", label = TRUE, pt.size = .1)

plot_grid(p1, p2)

new.cluster.ids <- c("Early_FGC","Late_FGC","DLC", "Early_FGC",
                       "LP", "Late_FGC","SC")

names(new.cluster.ids) <- levels(sample.all.qc.anchor)
sample.all.qc.anchor <- RenameIdents(sample.all.qc.anchor, new.cluster.ids)
sample.all.qc.anchor@meta.data$identi <- sample.all.qc.anchor@active.ident

DefaultAssay(sample.all.qc.anchor) <- "RNA"
anchored.markers <- FindAllMarkers(object = sample.all.qc.anchor, logfc.threshold=1.5, min.pct=0.6, only.pos = TRUE,test.use = "roc")
anchored.markers %>% group_by(cluster)
anchored_top10 <-anchored.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_diff)

DoHeatmap(object = subset(sample.all.qc.anchor,downsample = 500), features = anchored_top10$gene,group.bar = TRUE,
          size=4,angle=45,draw.lines=T)

# Marker gene plot
Early_FGC_marker<-c("NANOG","POU5F1","LIN28A","TFAP2C","NANOS3")
Late_FGC_marker<-c("DDX4","SIX1","TEX15","PIWIL2","PIWIL4")
LP_marker<-c("DCN", "COL1A1")
DLC_marker<-c("STAR","CYP17A1")
SC<-c("AMH","CITED1","SOX9")
plot_features <- c(LP_marker,DLC_marker,SC,
                   Early_FGC_marker,Late_FGC_marker)

DotPlot(sample.all.qc.anchor, features = plot_features)+
  RotatedAxis()

# get  DEG of FGCs

Idents(sample.all.qc.anchor) <- sample.all.qc.anchor@meta.data$identi

sample.all.qc.anchor@meta.data$celltype_now<-paste(sample.all.qc.anchor@meta.data$identi,sample.all.qc.anchor@meta.data$sex_type,sep = "_")
Idents(sample.all.qc.anchor) <- sample.all.qc.anchor@meta.data$celltype_now


Early_FGC_DEG <- FindMarkers(object = sample.all.qc.anchor, ident.1 = "Early_FGC_KS", ident.2 = "Early_FGC_Ctrl_M",
                             logfc.threshold = log(2), min.pct = 0.25,test.use ="roc",only.pos=F)%>%.[.$power>=0.4,]

Late_FGC_DEG <- FindMarkers(object = sample.all.qc.anchor, ident.1 = "Late_FGC_KS", ident.2 = "Late_FGC_Ctrl_M",
                            logfc.threshold = log(2), min.pct = 0.25,test.use ="roc",only.pos=F)%>%.[.$power>=0.4,]




 
