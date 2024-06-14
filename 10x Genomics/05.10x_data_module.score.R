## load seurat object
ser_anchor=readRDS("ser_anchor.rds")
Idents(ser_anchor)=factor(ser_anchor$celltype,levels = c("SPC","SC1","SC2"), ordered = TRUE)

## load gene set for module score
migration=read_excel("Genes_related_to_cell_migration2.29.xlsx",sheet = 2) 
adhesion=read_excel("Genes_related_to_cell_adhesion3.5.xlsx",sheet = 2)
mco=read_excel("Genes_related_to_microtubule_cytoskeleton_organization3.5.xlsx",sheet = 2)
cc=read_excel("Table_S3_Gene_list_for_GSEA_analysis.xlsx",sheet = 2)
qc.ser.nc.markers=read_xls("Sertoli.marker.in_Ctrl_M.xls")  

migration.gene=intersect(rownames(ser_anchor),migration$`Gene/product (bioentity_label)`) 
adhesion.gene=intersect(rownames(ser_anchor),adhesion$`Gene/product (bioentity_label)`) 
mco.gene=intersect(rownames(ser_anchor),mco$`Gene/product (bioentity_label)`) 
cc.gene=intersect(rownames(ser_anchor),union(cc$`G1/S`,cc$`G2/M`)) 
spcmarker.gene=unique(qc.ser.nc.markers[qc.ser.nc.markers$cluster=="SPC",]$gene) 
sc1marker.gene=unique(qc.ser.nc.markers[qc.ser.nc.markers$cluster=="SC1",]$gene) 
sc2marker.gene=unique(qc.ser.nc.markers[qc.ser.nc.markers$cluster=="SC2",]$gene) 

## calculate and get a Seurat object with module scores added to object meta data
ser_anchor <- AddModuleScore( 
  object = ser_anchor,
  features =list(spcmarker.gene),
  name = 'spcmarker_sore'
)
ser_anchor <- AddModuleScore( 
  object = ser_anchor,
  features =list(sc1marker.gene),
  name = 'sc1marker_sore'
)
ser_anchor <- AddModuleScore( 
  object = ser_anchor,
  features =list(sc2marker.gene),
  name = 'sc2marker_sore'
)
ser_anchor <- AddModuleScore(
  object = ser_anchor,
  features =list(ser.gene),
  name = 'ser_sore'
)
ser_anchor <- AddModuleScore(
  object = ser_anchor,
  features =list(migration.gene),
  name = 'migration_sore'
)
ser_anchor <- AddModuleScore(
  object = ser_anchor,
  features =list(adhesion.gene),
  name = 'adhesion_sore'
)
ser_anchor <- AddModuleScore(
  object = ser_anchor,
  features =list(mco.gene),
  name = 'mco_sore'
)
ser_anchor <- AddModuleScore(
  object = ser_anchor,
  features =list(cc.gene),
  name = 'cc_sore'
)
## Calculate p-values for differences between sex within each cell type, and perform Benjamini-Hochberg multiple testing correction.
score=colnames(ser_anchor@meta.data)[13:20]
module_result=data.frame()
for(i in score){
  test_result <- compare_means(as.formula(paste(i, "~ sex")), group.by = "celltype", 
                               data = ser_anchor@meta.data, method = "wilcox.test", 
                               p.adjust.method = "BH")
  test_result$BH =p.adjust(test_result$p,method = "BH")
  module_result=rbind(module_result,test_result)
}
module_result$p.adj.signif <- ifelse(module_result$BH < 0.001, "***",
                                     ifelse(module_result$BH < 0.01, "**",
                                            ifelse(module_result$BH < 0.05, "*", "ns")))
fwrite(module_result,"module_10x_wilcox_P.adj_BH.xls",sep = "\t",quote = F)

## Boxplot displaying module score results with significant markers labeled based on BH correction results
ser_anchor@meta.data$sex=factor(ser_anchor@meta.data$sex,levels = c("Ctrl_M","KS"),ordered = T )
pdf(file = "module.pdf",height = 5.4,width =7 )
ggboxplot(ser_anchor@meta.data,x="celltype",y="spcmarker_sore1",fill ="sex" ,outlier.shape = NA,palette = c("#00AFBB","#E7B800"))+
  stat_compare_means(aes(group =sex), label = "p.signif", method = "wilcox.test",label.y = 0.75,size = 6)+
  rotate_x_text(45)+
  ggtitle("SPC marker in Ctrl_M")+
  theme(legend.text = element_text(size = 12))+
  labs(y = "Module Score")
ggboxplot(ser_anchor@meta.data,x="celltype",y="sc1marker_sore1",fill ="sex" ,outlier.shape = NA,palette = c("#00AFBB","#E7B800"))+
  stat_compare_means(aes(group =sex), label = "p.signif", method = "wilcox.test",label.y = 2.4,size = 6)+
  rotate_x_text(45)+
  ggtitle("SC1 marker in Ctrl_M")+
  theme(legend.text = element_text(size = 12))+
  labs(y = "Module Score")
ggboxplot(ser_anchor@meta.data,x="celltype",y="sc2marker_sore1",fill ="sex" ,outlier.shape = NA,palette = c("#00AFBB","#E7B800"))+
  stat_compare_means(aes(group =sex), label = "p.signif", method = "wilcox.test",label.y = 0.7,size = 6)+
  rotate_x_text(45)+
  ggtitle("SC2 marker in Ctrl_M")+
  theme(legend.text = element_text(size = 12))+
  labs(y = "Module Score")
ggboxplot(ser_anchor@meta.data,x="celltype",y="ser_sore1",fill ="sex" ,outlier.shape = NA,palette = c("#00AFBB","#E7B800"))+
  stat_compare_means(aes(group =sex),label =  "p.signif", method = "wilcox.test",label.y = 0.6,size = 6)+
  rotate_x_text(45)+
  ggtitle("Sertoli marker")+
  theme(legend.text = element_text(size = 12))+
  labs(y = "Module Score")
ggboxplot(ser_anchor@meta.data,x="celltype",y="migration_sore1",fill ="sex" ,outlier.shape = NA,palette = c("#00AFBB","#E7B800"))+
  stat_compare_means(aes(group =sex),label =  "p.signif", method = "wilcox.test",label.y = 0.15,size = 6)+
  rotate_x_text(45)+
  ggtitle("Genes related to cell migration")+
  theme(legend.text = element_text(size = 12))+
  labs(y = "Module Score")
ggboxplot(ser_anchor@meta.data,x="celltype",y="adhesion_sore1",fill ="sex" ,outlier.shape = NA,palette = c("#00AFBB","#E7B800"))+
  stat_compare_means(aes(group =sex),label =  "p.signif", method = "wilcox.test",label.y = 0.1,size = 6)+
  rotate_x_text(45)+
  ggtitle("Genes related to cell adhesion")+
  theme(legend.text = element_text(size = 12))+
  labs(y = "Module Score")
ggboxplot(ser_anchor@meta.data,x="celltype",y="mco_sore1",fill ="sex" ,outlier.shape = NA,palette = c("#00AFBB","#E7B800"))+
  stat_compare_means(aes(group =sex),label =  "p.signif", method = "wilcox.test",label.y = 0.25,size = 6)+
  rotate_x_text(45)+
  ggtitle("Genes related to microtubule cytoskeleton organization")+
  theme(legend.text = element_text(size = 12))+
  labs(y = "Module Score")
ggboxplot(ser_anchor@meta.data,x="celltype",y="cc_sore1",fill ="sex" ,outlier.shape = NA,palette = c("#00AFBB","#E7B800"))+
  stat_compare_means(aes(group =sex),label =  "p.signif", method = "wilcox.test",label.y = 0.85,size = 6)+
  rotate_x_text(45)+
  ggtitle("Genes related to cell cycle")+
  theme(legend.text = element_text(size = 12))+
  labs(y = "Module Score")
dev.off()
