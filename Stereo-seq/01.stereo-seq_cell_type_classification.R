options(stringsAsFactors = FALSE)
library(readxl)
library(data.table)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(patchwork)
library(Seurat)
library(Matrix)
library(ggpubr)
library(ggrepel)
library(openxlsx)
library(future)
## read stereo-seq data matrix
ST.data.TNC16 <- read.table("TNC16_D03058F1_hg19/D03058F1.tissue.gem", sep = "\t", header = T, check.names = F)# read data matrix
ST.data.TNC18 <- read.table("TNC18_D03058E6_hg19/D03058E6.tissue.gem", sep = "\t", header = T, check.names = F)# read data matrix
ST.data.TKS5 <- read.table("TKS5_SS200000513TL_C5_hg19/SS200000513TL_C5.tissue.gem", sep = "\t", header = T, check.names = F)# read data matrix
ST.data.TKS6 <- read.table("TKS6_SS200000513TL_C6_hg19/SS200000513TL_C6.tissue.gem", sep = "\t", header = T, check.names = F)# read data matrix
## define a function to Create Stereo-seq object based on DoubletFinder package documentation.
Creat.ST.object <- function(ST.data, bin){
  ST.data.bin40 <- ST.data
  ST.data.bin40$x1 <- trunc(ST.data.bin40$x/bin) * bin
  ST.data.bin40$y1 <- trunc(ST.data.bin40$y/bin) * bin
  ST.data.bin40$cellID <- paste(ST.data.bin40$x1, "_", ST.data.bin40$y1, sep = "")
  ST.data.bin40 <- aggregate(ST.data.bin40$MIDCount, by = list(ST.data.bin40$cellID, ST.data.bin40$geneID), sum)
  colnames(ST.data.bin40) <- c("cellID", "geneID", "MIDCounts")
  ST.data.bin40$cellInx <- match(ST.data.bin40$cellID, unique(ST.data.bin40$cellID))
  ST.data.bin40$cellInx <- match(ST.data.bin40$cellID, unique(ST.data.bin40$cellID))
  ST.data.bin40$geneInx <- match(ST.data.bin40$geneID, unique(ST.data.bin40$geneID))
  mat <- sparseMatrix(i = ST.data.bin40$geneInx, j = ST.data.bin40$cellInx, x = ST.data.bin40$MIDCounts, 
                      dimnames = list(unique(ST.data.bin40$geneID), unique(ST.data.bin40$cellID)))
  ST.bin40.coord.df <- data.frame(cellname = colnames(mat))
  rownames(ST.bin40.coord.df) <- ST.bin40.coord.df$cellname
  ST.bin40.coord.df <- separate(ST.bin40.coord.df, col = cellname, sep = "_", into = c("x", "y"))
  ST.bin40 <- CreateSeuratObject(mat, project =paste0("ST.bin",bin), assay = "Spatial")
  ST.bin40$slice <- 1
  ST.bin40$region <- paste0("ST.bin",bin)
  colnames(ST.bin40.coord.df) <- c("imagerow", "imagecol")
  ST.bin40.coord.df$imagerow <- as.numeric(ST.bin40.coord.df$imagerow)
  ST.bin40.coord.df$imagecol <- as.numeric(ST.bin40.coord.df$imagecol)
  ST.bin40@images$ST_bin40 <- new(Class = "SlideSeq", assay = "spatial", key = "image_", coordinates = ST.bin40.coord.df)
  ## SCT transform
  ST.bin40 <- SCTransform(ST.bin40, assay = "Spatial", verbose = FALSE)
  ST.bin40 <- RunPCA(ST.bin40, verbose = F)
  ElbowPlot(ST.bin40, ndims = 30)
  ST.bin40 <- FindNeighbors(ST.bin40, reduction = "pca", dims = 1:20)
  ST.bin40 <- FindClusters(ST.bin40, verbose = T, resolution = 0.5)
  ST.bin40 <- RunUMAP(ST.bin40, reduction = "pca", dims = 1:20)
  return(ST.bin40)
}
## load 10x scRNA-seq data and SCT conversion
sample.finally=readRDS("sample.finally.rds")
sample.finally.sct <- SCTransform(sample.finally,ncells = 3000,verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30) 
Idents(sample.finally.sct)=sample.finally.sct$celltype

#######################Assign the cell type for each spatial spot with different bin
for(b in c(10,20,30,40,50,100)){
  setwd(paste0("./bin",b))
  for(i in c("TNC16","TNC18","TKS5","TKS6")){
    gem.name=paste0("ST.data.",i)
    gem.data <- get(gem.name)
    ## Creat the seurat object for each data on each slide
    gem=Creat.ST.object(gem.data, b)
    ## integrate scRNAs-seq and ST.seurat.object
    anchors<- FindTransferAnchors(reference = sample.finally.sct, query = gem,normalization.method = "SCT")
    predictions <- TransferData(anchorset =anchors, refdata = sample.finally.sct$celltype, prediction.assay = F,
                                weight.reduction = gem[["pca"]], dims = 1:30)
    predic_mtrx <- as.data.frame(predictions)
    gem$Celltype <- predic_mtrx$predicted.id
    predictions.assay <- TransferData(anchorset = anchors, refdata = sample.finally.sct$celltype, 
                                      prediction.assay = TRUE, weight.reduction = gem[["pca"]], dims = 1:30)
    ## Assign the cell type for each spatial spot
    gem[["predictions"]] <- predictions.assay
    DefaultAssay(gem) <- "predictions"
    celltype.score <- GetAssayData(gem, slot = "data")
    celltype.score <-as.data.frame(celltype.score)
    celltype.max <- apply(celltype.score, 2, function(x) rownames(celltype.score)[which.max(x)])
    celltype <- data.frame(cell_type = celltype.max, 
                           cell_name = names(celltype.max), 
                           GetTissueCoordinates(gem), 
                           spatial_clusters = gem$seurat_clusters,
                           nFeature_RNA = gem$nFeature_Spatial, 
                           nCount_RNA = gem$nCount_Spatial)
    score=as.data.frame(t(celltype.score[15,]))
    gem@meta.data$score=score$max
    celltype.seurat <- CreateSeuratObject(counts = celltype.score[-15,], project = "celltype", min.cells = 0, min.features = 0,meta.data = celltype)
    table(celltype.seurat$cell_type)
    celltype.seurat$cell_type=factor(celltype.seurat$cell_type,levels =c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","SPC","SC1","SC2","PMC","EC","MC","NK"),ordered = T )
    p=DoHeatmap(celltype.seurat,features = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","SPC","SC1","SC2","PMC","EC","MC","NK"),assay = "RNA",slot = "counts" ,group.by ="cell_type")+
      scale_fill_gradientn(colors = c("white","firebrick3"))
    ggsave(paste0(i,".bin",b,".sp.ct.pdf"), p, width = 25, height = 15, device = "pdf")
    celltype.score$celltype=rownames(celltype.score)
    celltype.score=celltype.score[, c(ncol(celltype.score), 1:(ncol(celltype.score)-1))]
    write.table(celltype.score, file = paste0(i,".bin",b,".celltype.score.xls"), quote=F, sep="\t", col.names=T, row.names=F)
    gem@meta.data$celltype=celltype.max
    gem$celltype=factor(gem$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    Idents(gem)=gem@meta.data$celltype
    score.sample=paste0(i,".bin",b,".celltype.score")
    SCT.object=paste0("SCT.",i,".bin",b)
    assign(SCT.object, gem)
    assign(score.sample, celltype.score)
    saveRDS(gem,file = paste0("SCT.",i,".bin",b,".rds"))
  }
}
################################select suitable bin
## Visualize the distribution of nCount and nFeature for each bin
for(b in c(10,20,30,40,50,100)){
  for(i in c("TNC16","TNC18","TKS5","TKS6")){
    SCT.object=paste0("SCT.",i,".bin",b)
    seurat <- get(SCT.object)
    plot1 <- VlnPlot(seurat, features = "nCount_Spatial", pt.size = 0,group.by = "orig.ident") + NoLegend()
    plot2 <- VlnPlot(seurat, features = "nFeature_Spatial", pt.size = 0,group.by = "orig.ident") + NoLegend()
    plot_grid(plot1, plot2,ncol = 2)
    plot3 <- SpatialFeaturePlot(seurat, features = "nCount_Spatial") + theme(legend.position = "right")
    plot4 <- SpatialFeaturePlot(seurat, features = "nFeature_Spatial") + theme(legend.position = "right")
    plot_grid(plot3, plot4,ncol = 2)
  }
}
## Calculate the number of each cell type per bin for each sample                         
nct=data.frame()
for(b in c(10,20,30,40,50,100)){
  for(i in c("TNC16","TNC18","TKS5","TKS6")){
    SCT.object=paste0("SCT.",i,".bin",b)
    SCT.data <- get(SCT.object)
    SCT.data$celltype=factor(SCT.data$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    print(SCT.object) 
    print(table(SCT.data$celltype))
    ct=as.data.frame(table(SCT.data$celltype))
    ct$sample=rep(i,nrow(ct))
    ct$bin=rep(b,nrow(ct))
    nct=rbind(nct,ct)
  }
}
write.table(nct, file = "celltype_count.xls", quote=F, sep="\t", col.names=T, row.names=F)

## Overall comparison shows that 40 is a suitable bin
