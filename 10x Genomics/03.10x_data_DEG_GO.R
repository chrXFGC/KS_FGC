###############03.We used FindMarkers function to calculate DEG and used clusterProfiler package to perform GO enrichment analysis
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
######################### calculate differentially expressed genes (DEGs) of each cell type between male controls and KS
sample.finally@meta.data$CT_Batch=NA
sample.finally@meta.data$CT_Batch <- paste(sample.finally@meta.data$celltype, sample.finally@meta.data$sex, sep = "_")
table(sample.finally@meta.data$CT_Batch)
sample.finally@meta.data$CT_Batch <-factor(sample.finally@meta.data$CT_Batch,levels =c("Early_FGC_KS","Early_FGC_Ctrl_M","Late_FGC_KS","Late_FGC_Ctrl_M","LPC1_KS","LPC1_Ctrl_M","LPC2_KS","LPC2_Ctrl_M","LP1_KS","LP1_Ctrl_M","LP2_KS","LP2_Ctrl_M","DLC_KS","DLC_Ctrl_M","PMC_KS","PMC_Ctrl_M","SPC_KS","SPC_Ctrl_M","SC1_KS","SC1_Ctrl_M","SC2_KS","SC2_Ctrl_M","EC_KS","EC_Ctrl_M","MC_KS","MC_Ctrl_M","NK_KS","NK_Ctrl_M"),ordered = T)
     
DEG.result<- data.frame()
Idents(sample.finally) <- sample.finally@meta.data$CT_Batch
levels(sample.finally)
for(i in c("Early_FGC","Late_FGC","LPC1","LPC2","LP1","LP2","DLC","PMC","SPC","SC1","SC2","EC","MC","NK")) {
  i1=paste(i, "KS", sep = "_")
  i2=paste(i, "NC", sep = "_")
  qcmarker<- FindMarkers(sample.finally, ident.1 = i1, ident.2 = i2, only.pos = F, min.pct = 0.3, logFC.thresh = 0.25)
  qcmarker=qcmarker[qcmarker$p_val<0.05,] 
  qcmarker$gene=row.names(qcmarker)
  qcmarker$state=paste("DEG",i, "KS_vs_NC", sep = "_")
  qcmarker$express=ifelse(qcmarker$avg_logFC >0, "KS_high","KS_low")
  DEG.result=rbind(DEG.result,qcmarker)
}
write.table(DEG.result, "DEG.result.qc.xls", quote=F, sep="\t", col.names=T, row.names=F)

########################## GO enrichment analysis on differentially expressed genes
unique_states <- unique(DEG.result$state)
GO.result<- data.frame()
go_results_list<- list()
for(i in unique_states) {
  deglist_high=DEG.result$gene[DEG.result$state == i&DEG.result$express == "KS_high"]
  erich.go.BP = enrichGO(gene =deglist_high,
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",
                         pvalueCutoff = 0.05)
  state <- sub("^[^_]+_", "", i)
  go_results_list[[paste("GO","KS_high", state, sep = "_")]] <- erich.go.BP
  erich.go.BP@result$state=rep(paste("KS_high",state,sep = "_"),nrow(erich.go.BP@result))
  GO.result=rbind(GO.result,erich.go.BP@result)
  
  deglist_low=DEG.result$gene[DEG.result$state == i&DEG.result$express == "KS_low"]
  if (length(deglist_low)> 0) {
    erich.go.BP2 = enrichGO(gene =deglist_low,
                            OrgDb = org.Hs.eg.db,
                            keyType = "SYMBOL",
                            ont = "BP",
                            pvalueCutoff = 0.05)
    state <- sub("^[^_]+_", "", i)
    go_results_list[[paste("GO","KS_low", state, sep = "_")]] <- erich.go.BP2
    erich.go.BP2@result$state=rep(paste("KS_low",state,sep = "_"),nrow(erich.go.BP2@result))
    GO.result=rbind(GO.result,erich.go.BP2@result)
  }
  else {
    next
  }
}
write.table(GO.result, "GO.result.qc.xls", quote=F, sep="\t", col.names=T, row.names=F)

## Save the results into separate sheets in xlsx format
goresultexcel <- createWorkbook()
for (i in 1:length(go_results_list)) {
  resultn <- go_results_list[[i]]
  sheet_name <- names(go_results_list)[i]
  sheet_name <- sub("^[^_]+_", "", sheet_name)
  addWorksheet(goresultexcel, sheet_name)
  writeData(goresultexcel, sheet = sheet_name, x = resultn@result)
}
saveWorkbook(goresultexcel, file = "go_result_list_qc.xlsx")

########################## GO enrichment analysis on marker genes found in sertoli cells of male controls
qc.ser.nc.markers=read_xls("Sertoli.marker.in_Ctrl_M.xls")
unique_states.ser.nc <- unique(qc.ser.nc.markers$cluster)
ser.nc.GO.result<- data.frame()
ser.nc.go_results_list<- list()

for(i in unique_states.ser.nc) {
  deglist=qc.ser.nc.markers$gene[qc.ser.nc.markers$cluster == i]
  erich.go.BP = enrichGO(gene =deglist,
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",
                         pvalueCutoff = 0.05)
  #state <- sub("^[^_]+_", "", i)
  ser.nc.go_results_list[[paste("GO", i, sep = "_")]] <- erich.go.BP
  erich.go.BP@result$state=rep(i,nrow(erich.go.BP@result))
  ser.nc.GO.result=rbind(ser.nc.GO.result,erich.go.BP@result)
}
write.table(ser.nc.GO.result, "sertoli.nc.GO.result.xls", quote=F, sep="\t", col.names=T, row.names=F)

## Save the results into separate sheets in xlsx format
goresultexcel <- createWorkbook()
for (i in 1:length(ser.nc.go_results_list)) {
  resultn <- ser.nc.go_results_list[[i]]
  sheet_name <- names(ser.nc.go_results_list)[i]
  #sheet_name <- sub("^[^_]+_", "", sheet_name)
  addWorksheet(goresultexcel, sheet_name)
  writeData(goresultexcel, sheet = sheet_name, x = resultn@result)
}
saveWorkbook(goresultexcel, file = "sertoli.nc.go_result_list.xlsx")

