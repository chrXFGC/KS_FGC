
ST.data.TNC16 <- read.table("~/analysis/KS/steroseq/01.gem/TNC16_D03058F1_hg19/D03058F1.tissue.gem", sep = "\t", header = T, check.names = F)# read data matrix
ST.data.TNC18 <- read.table("~/analysis/KS/steroseq/01.gem/TNC18_D03058E6_hg19/D03058E6.tissue.gem", sep = "\t", header = T, check.names = F)# read data matrix
ST.data.TKS5 <- read.table("~/analysis/KS/steroseq/01.gem/TKS5_SS200000513TL_C5_hg19/SS200000513TL_C5.tissue.gem", sep = "\t", header = T, check.names = F)# read data matrix
ST.data.TKS6 <- read.table("~/analysis/KS/steroseq/01.gem/TKS6_SS200000513TL_C6_hg19/SS200000513TL_C6.tissue.gem", sep = "\t", header = T, check.names = F)# read data matrix

Creat.ST.object <- function(ST.data, bin,nslice,region){
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
  ST.bin40$slice <- nslice
  ST.bin40$region <- region
  colnames(ST.bin40.coord.df) <- c("imagerow", "imagecol")
  ST.bin40.coord.df$imagerow <- as.numeric(ST.bin40.coord.df$imagerow)
  ST.bin40.coord.df$imagecol <- as.numeric(ST.bin40.coord.df$imagecol)
  ST.bin40@images$ST_bin40 <- new(Class = "SlideSeq", assay = "spatial", key = "image_", coordinates = ST.bin40.coord.df)
  return(ST.bin40)
}
SCT.TKS5.merge=Creat.ST.object(ST.data = ST.data.TKS5,bin = 40,nslice = 1,region = "ST.TKS5")
SCT.TKS6.merge=Creat.ST.object(ST.data = ST.data.TKS6,bin = 40,nslice = 2,region = "ST.TKS6")
SCT.TNC16.merge=Creat.ST.object(ST.data = ST.data.TNC16,bin = 40,nslice = 3,region = "ST.TNC16")
SCT.TNC18.merge=Creat.ST.object(ST.data = ST.data.TNC18,bin = 40,nslice = 4,region = "ST.TNC18")

for(i in c("TNC16","TNC18","TKS5","TKS6")){
  SCT.object=paste0("SCT.",i,".bin40")
  seurat <- get(SCT.object)
  t1=as.data.frame(seurat@assays[["predictions"]]@data)
  t2=as.data.frame(t(t1[15,]))
  seurat@meta.data$max=t2$max
  meta=seurat@meta.data
  meta$cell=row.names(meta)
  meta=meta[,c(10,11,12,13)]
  merge.object==paste0("SCT.",i,".merge")
  seurat2 <- get(merge.object)
  seurat2@meta.data$cell=row.names(seurat2@meta.data)
  seurat2@meta.data$Celltype=NA
  seurat2@meta.data$celltype=NA
  seurat2@meta.data$max=NA
  for(i in 1:nrow(meta18)){
    seurat2@meta.data[which(seurat2@meta.data$cell == meta$cell[i]),'Celltype'] <- as.character(meta$Celltype[i])
    seurat2@meta.data[which(seurat2@meta.data$cell == meta$cell[i]),'celltype'] <- as.character(meta$celltype[i])
    seurat2@meta.data[which(seurat2@meta.data$cell == meta$cell[i]),'max'] <- as.character(meta$max[i])
  }
  colnames(seurat2)=paste0(i,"_", colnames(seurat2))
  SCT.data$celltype=factor(SCT.data$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
  Idents(SCT.data)=factor(SCT.data$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
  assign(merge.object,seurat2)
}

ST.merge=merge(SCT.TKS5.merge,c(SCT.TKS6.merge,SCT.TNC16.merge,SCT.TNC18.merge))
ST.merge <- SCTransform(ST.merge, assay = "Spatial", verbose = FALSE)

c16=ST.merge[,row.names(ST.merge@meta.data[ST.merge@meta.data$region=="ST.TNC16",])]
c18=ST.merge[,row.names(ST.merge@meta.data[ST.merge@meta.data$region=="ST.TNC18",])]
c5=ST.merge[,row.names(ST.merge@meta.data[ST.merge@meta.data$region=="ST.TKS5",])]
c6=ST.merge[,row.names(ST.merge@meta.data[ST.merge@meta.data$region=="ST.TKS6",])]

for(i in c("16","18","6","5")){ 
  sct=get(paste0("c",i))
  Idents(sct)=factor(sct$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
  if(i %in% c("16")){
    SpatialDimPlot(sct, crop = TRUE,cols = color.use2)
    c161=sct[,rownames(sct@images[["ST_bin40.3"]]@coordinates[sct@images[["ST_bin40.3"]]@coordinates$x<= 6960 & sct@images[["ST_bin40.3"]]@coordinates$y >= 12920,])]
    Idents(c161)=factor(c161$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    c161@meta.data$ident=rep("Ctrl_M8_1",ncol(c161))
    SpatialDimPlot(c161, crop = TRUE,cols = color.use2,stroke = 0)
    c162=sct[,rownames(sct@images[["ST_bin40.3"]]@coordinates[sct@images[["ST_bin40.3"]]@coordinates$x>= 12280 & sct@images[["ST_bin40.3"]]@coordinates$y <= 14040,])]
    Idents(c162)=factor(c162$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    SpatialDimPlot(c162, crop = TRUE,cols = color.use2)
    c162@meta.data$ident=rep("Ctrl_M8_2",ncol(c162))
  }
  if(i %in% c("18")){
    SpatialDimPlot(sct, crop = TRUE,cols = color.use2)
    c181=sct[,rownames(sct@images[["ST_bin40.4"]]@coordinates[sct@images[["ST_bin40.4"]]@coordinates$x<= 8720 & sct@images[["ST_bin40.4"]]@coordinates$y >= 12800,])]
    Idents(c181)=factor(c181$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    SpatialDimPlot(c181, crop = TRUE,cols = color.use2)
    c181@meta.data$ident=rep("Ctrl_M9_1",ncol(c181))
    c182=sct[,rownames(sct@images[["ST_bin40.4"]]@coordinates[sct@images[["ST_bin40.4"]]@coordinates$x>= 12480 & sct@images[["ST_bin40.4"]]@coordinates$y <= 11120,])]
    Idents(c182)=factor(c182$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    SpatialDimPlot(c182, crop = TRUE,cols = color.use2)
    c182@meta.data$ident=rep("Ctrl_M9_2",ncol(c182))
  }
  if(i %in% c("6")){
    SpatialDimPlot(sct, crop = TRUE,cols = color.use2)
    c61=sct[,rownames(sct@images[["ST_bin40.2"]]@coordinates[sct@images[["ST_bin40.2"]]@coordinates$x<= 12080 & sct@images[["ST_bin40.2"]]@coordinates$y <= 16160,])]
    Idents(c61)=factor(c61$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    SpatialDimPlot(c61, crop = TRUE,cols = color.use2)
    c61@meta.data$ident=rep("KS13_2",ncol(c61))
    c62=sct[,rownames(sct@images[["ST_bin40.2"]]@coordinates[sct@images[["ST_bin40.2"]]@coordinates$x>= 18560 & sct@images[["ST_bin40.2"]]@coordinates$y >= 11760,])]
    Idents(c62)=factor(c62$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    SpatialDimPlot(c62, crop = TRUE,cols = color.use2)
    c62@meta.data$ident=rep("KS13_2",ncol(c62))
  }
  if(i %in% c("5")){
    SpatialDimPlot(sct, crop = TRUE,cols = color.use2)
    c5_1=sct[,rownames(sct@images[["ST_bin40"]]@coordinates[sct@images[["ST_bin40"]]@coordinates$x<= 14080 & sct@images[["ST_bin40"]]@coordinates$y <= 11920,])]
    Idents(c5_1)=factor(c5_1$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    SpatialDimPlot(c5_1, crop = TRUE,cols = color.use2)
    c51@meta.data$ident=rep("KS12_1",ncol(c51))
    c5_2=sct[,rownames(sct@images[["ST_bin40"]]@coordinates[sct@images[["ST_bin40"]]@coordinates$x<= 10920 & sct@images[["ST_bin40"]]@coordinates$y >= 17760,])]
    Idents(c5_2)=factor(c5_2$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    SpatialDimPlot(c5_2, crop = TRUE,cols = color.use2)
    c52@meta.data$ident=rep("KS12_2",ncol(c52))
    c5_3=sct[,rownames(sct@images[["ST_bin40"]]@coordinates[sct@images[["ST_bin40"]]@coordinates$x>= 18280 & sct@images[["ST_bin40"]]@coordinates$y >= 17600,])]
    Idents(c5_3)=factor(c5_3$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
    SpatialDimPlot(c5_3, crop = TRUE,cols = color.use2)
    c53@meta.data$ident=rep("KS12_3",ncol(c53))
  }
}

c5_1.c=c5_1[,rownames(c5_1@images[["ST_bin40"]]@coordinates[c5_1@images[["ST_bin40"]]@coordinates$y <= 9320 &
                                                              c5_1@images[["ST_bin40"]]@coordinates$y >= 8500 &
                                                              c5_1@images[["ST_bin40"]]@coordinates$x <= 12160 &
                                                              c5_1@images[["ST_bin40"]]@coordinates$x >= 11320,])]
c51.delet=row.names(c5_1.c@meta.data[c5_1.c@meta.data$max<0.3,]) #136
c51=c5_1[,setdiff(colnames(c5_1),c51.delet)]
Idents(c51)=factor(c51$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
SpatialDimPlot(c51, crop = TRUE,cols = color.use2)+
  geom_vline(xintercept = seq(7000, 12000, by = 500), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(8500, 14500, by = 500), color = "black", linetype = "dashed")
c5_2.c=c5_2[,rownames(c5_2@images[["ST_bin40"]]@coordinates[c5_2@images[["ST_bin40"]]@coordinates$x <= 8480 &
                                                              c5_2@images[["ST_bin40"]]@coordinates$x >= 7480 &
                                                              c5_2@images[["ST_bin40"]]@coordinates$y <= 21420 &
                                                              c5_2@images[["ST_bin40"]]@coordinates$y >= 20080,])]
c52.delet=row.names(c5_2.c@meta.data[c5_2.c@meta.data$max<0.3,]) #136
c52=c5_2[,setdiff(colnames(c5_2),c52.delet)]
Idents(c52)=factor(c52$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
SpatialDimPlot(c52, crop = TRUE,cols = color.use2)+
  geom_vline(xintercept = seq(17500, 22500, by = 500), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(5500, 11000, by = 500), color = "black", linetype = "dashed")
c5_3.c=c5_3[,rownames(c5_3@images[["ST_bin40"]]@coordinates[c5_3@images[["ST_bin40"]]@coordinates$y <= 20180 &
                                                              c5_3@images[["ST_bin40"]]@coordinates$y >= 19220 &
                                                              c5_3@images[["ST_bin40"]]@coordinates$x <= 21900 &
                                                              c5_3@images[["ST_bin40"]]@coordinates$x >= 20480,])]
c53.delet=row.names(c5_3.c@meta.data[c5_3.c@meta.data$max<0.3,]) #136
c53=c5_3[,setdiff(colnames(c5_3),c53.delet)]
Idents(c53)=factor(c53$celltype,levels = c("Early-FGC","Late-FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK"), ordered = TRUE)
SpatialDimPlot(c53, crop = TRUE,cols = color.use2)+
  geom_vline(xintercept = seq(17600, 22500, by = 500), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(18000, 23000, by = 500), color = "black", linetype = "dashed")

saveRDS("c161",file = "Ctrl_M8_1.rds")
saveRDS("c181",file = "Ctrl_M9_1.rds")
saveRDS("c51",file = "KS12_1.rds")
saveRDS("c61",file = "KS13_1.rds")
saveRDS("c162",file = "Ctrl_M8_2.rds")
saveRDS("c182",file = "Ctrl_M9_2.rds")
saveRDS("c52",file = "KS12_2.rds")
saveRDS("c62",file = "KS13_2.rds")
saveRDS("c53",file = "KS12_3.rds")