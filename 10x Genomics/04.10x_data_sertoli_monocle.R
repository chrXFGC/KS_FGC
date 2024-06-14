##############################monocle
library(monocle)
####################sertoli
## load up sertoli data into Monocle's main class, CellDataSet
ser_anchor=readRDS("ser_anchor.rds")
expr_matrix.ser <- as(as.matrix(ser_anchor@assays$RNA@counts), 'sparseMatrix')
dim(expr_matrix.ser)
p_data.ser <- ser_anchor@meta.data
f_data.ser <- data.frame(gene_short_name = row.names(expr_matrix.ser),row.names = row.names(expr_matrix.ser)) 

pd.ser <- new('AnnotatedDataFrame', data = p_data.ser)
fd.ser <- new('AnnotatedDataFrame', data = f_data.ser)
cds.ser <- newCellDataSet(expr_matrix.ser,
                          phenoData = pd.ser,
                          featureData = fd.ser,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial())
## Estimate size factors and dispersions
cds.ser <- estimateSizeFactors(cds.ser)
cds.ser <- estimateDispersions(cds.ser)
## Filtering low-quality cells
cds.ser <- detectGenes(cds.ser, min_expr = 0.1) 
print(head(fData(cds.ser)))
expressed_genes.ser <- row.names(subset(fData(cds.ser),
                                        num_cells_expressed >= 10)) 
print(head(expressed_genes.ser))
## choose genes that define a cell's progress
ser.markers=read_xls("Sertoli.marker.xls")
cds.ser <- setOrderingFilter(cds.ser, ser.markers$gene)
plot_ordering_genes(cds.ser)
## reduce data dimensionality
cds.ser <- reduceDimension(cds.ser, max_components = 2,method = 'DDRTree', verbose=T, sigma = 0.1)
## order cells along the trajectory
cds.ser <- orderCells(cds.ser, reverse = T)
## visualize the trajectory in the reduced dimensional space
plot_cell_trajectory(cds.ser, color_by = "State",show_branch_points=F)+ facet_wrap("~sex", nrow = 1)
plot_cell_trajectory(cds.ser, color_by = "Pseudotime",show_branch_points=F)+ facet_wrap("~batch", nrow = 1)
plot_cell_trajectory(cds.ser, color_by = "celltype",show_branch_points=F)+ facet_wrap("~batch", nrow = 1)
plot_cell_trajectory(cds.ser, shape_by = "celltype",show_branch_points=F, color_by = "Pseudotime")+ facet_wrap("~batch", nrow = 1)
plot_cell_trajectory(cds.ser, color_by = "Pseudotime",show_branch_points=F)
plot_cell_trajectory(cds.ser, shape_by = "Pseudotime",show_branch_points=F, color_by = "celltype")+ facet_wrap("~batch", nrow = 1)

plot_a<-plot_cell_trajectory(cds.ser, color_by ="celltype" ,show_backbone =  T,
                             show_tree = T,show_branch_points=F,theta=175)
##############################################plot
point.data <- ggplot_build(plot_a)[["plot"]][["data"]]
ser_color=c("#00BFFF","#87CEFA","#4b74b2")
names(ser_color)=c("SPC","SC1","SC2")
point.data$celltype <- factor(point.data$celltype, levels = c("SPC","SC1","SC2"), ordered = T)

######plot with alpha setting
ggscatter(point.data,x="data_dim_1",y="data_dim_2",color = "celltype",alpha=0.3,
              xlab = "Component 1",ylab = "Component 2",
              size = "Pseudotime",facet.by = "sex")+ #+ facet_wrap("~Sex", nrow = 1)+
  scale_color_manual(breaks = c( "SPC","SC1","SC2"), values=ser_color) +
  theme(legend.position = "right")+ scale_size_continuous(range = c(0.1, 3))+
  guides(color=guide_legend(""))
#######plot without alpha setting
ggscatter(point.data,x="data_dim_1",y="data_dim_2",color =  "celltype",
              xlab = "Component 1",ylab = "Component 2",
              size = "Pseudotime",facet.by = "sex")+ #+ facet_wrap("~Sex", nrow = 1)+
  scale_color_manual(breaks = c( "SPC","SC1","SC2"), values=ser_color) +
  theme(legend.position = "right")+ scale_size_continuous(range = c(0.1, 3))+
  guides(color=guide_legend(""))

######################umap plot
umap_result_coord <- as.data.frame(ser_anchor@reductions$umap@cell.embeddings)
umap_result_coord$Cluster <- ser_anchor@meta.data[rownames(umap_result_coord),"celltype"]
umap_result_coord$sample_name=rownames(umap_result_coord)

plot_df_corrd<-merge(point.data,umap_result_coord,by="sample_name")
ggscatter(plot_df_corrd,x="UMAP_1",y="UMAP_2",color =  "celltype",
          xlab = "UMAP_1",ylab = "UMAP_2",size = 1.2,
          facet.by = "sex",ncol = 1)+ #+ facet_wrap("~Sex", nrow = 1)+ #size = "Pseudotime",
  scale_color_manual(breaks = c( "SPC","SC1","SC2"), values=ser_color) +
  theme(legend.position = "right", 
        axis.ticks = element_blank())+ scale_size_continuous(range = c(0.1, 3))+
  guides(color=guide_legend(""))

ggscatter(plot_df_corrd, x = "UMAP_1", y = "UMAP_2", color = "Pseudotime",
          size = 1.2, facet.by = "sex",ncol = 1) +
  scale_color_gradient(low = "blue", high = "red", name = "Pseudotime") +  
  theme(legend.position = "right", 
        axis.ticks = element_blank()) 
