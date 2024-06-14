library(wesanderson)
library(plot3D)
library(readxl)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

xist_data=read_excel("XIST_scope_FINAL.xlsx",sheet = 1)
## Calculate the mean fluorescence intensity(MFI)
xist_data$av_l=xist_data$weighted_sum/xist_data$counts
## Calculate the despersion
xist_data$d=xist_data$area/xist_data$counts
xist_data$av_l[is.na(xist_data$av_l)] <- 0
xist_data$d[is.na(xist_data$d)] <- 0
## Calculate the mean fluorescence intensity of surrounding somatic cells on the same slide
average_soma=xist_data %>%
  filter(GROUP == "Soma") %>%  
  group_by(SAMPLE, MARKER) %>%  
  summarize(avl_mean = mean(av_l)) 
xist_data$avl_mean=NA
for (i in 1:nrow(average_soma)) {
  xist_data$avl_mean <- ifelse(xist_data$SAMPLE == average_soma$SAMPLE[i] & xist_data$MARKER == average_soma$MARKER[i],
                          average_soma$avl_mean[i], xist_data$avl_mean)
}
## Calculate the Relative MFI
xist_data$r_l <- xist_data$av_l/xist_data$avl_mean
## Delete the results for somatic cells
xist_data2=xist_data[xist_data$GROUP != "Soma",]
## Extract the sequence number for each cell, Area, Relative MFI, and Dispersion
data=data.frame(name=xist_data2$num,counts=xist_data2$counts,rl=xist_data2$r_l,d=xist_data2$d)
row.names(data)=data$name
data=data[,-1]
## scale data and perform k-means clustering with k=4
data <- scale(data)  
k <- 4  
kmeans_result <- kmeans(data, centers = k)
cluster_labels <- kmeans_result$cluster
result_with_cluster <- cbind(kmeans_result$layout, cluster = factor(cluster_labels))
result=as.data.frame(result_with_cluster)
result$name=rownames(result)
kmeans_result_with_cluster=cbind(result,xist_data2)
kmeans_result_with_cluster=kmeans_result_with_cluster[,-2]
## Calculate the number and proportion of each cluster
proportion <- kmeans_result_with_cluster %>%
  group_by(GROUP,MARKER,cluster) %>%
  tally() %>%
  transmute(GROUP,MARKER,cluster,n,lineage_sum=sum(n),pct = 100*n / sum(n))
proportion$cluster <- factor(proportion$cluster)
## Define and rename clusters based on the clustering results and features
proportion$state=ifelse(proportion$cluster==1,"no XIST",ifelse(proportion$cluster==2,"Nebular",ifelse(proportion$cluster==3,"Dispersed",ifelse(proportion$cluster==4,"Confined",""))))
proportion$state=factor(proportion$state,levels = c("no XIST", "Dispersed", "Confined", "Nebular"),ordered = T)
## Visualize the cell proportions of each state for each sex cell type facet by sex
proportion$MARKER=factor(proportion$MARKER,levels = c("OCT4","DDX4"),ordered = T)
proportion$GROUP=factor(proportion$GROUP,levels = c("KS","CTRL_F"),ordered = T) 
ggbarplot(proportion, x="MARKER", y="pct", fill = "state",color ="state",palette =c("#CDC1C5", "#96CDCD","#DE74B0", "#FFD39B"),
          xlab="Group",ylab="Percent",legend="right")+
  ggtitle("GROUP")+
  rotate_x_text(45)+guides(fill=guide_legend(title="GROUP"))+
  facet_wrap(~GROUP, nrow = 1)
## scatter plot for Area and Dispersion of each cell
ggplot(kmeans_result_with_cluster, aes(x = counts, y = d)) +
  geom_point(aes(size = r_l, color = state), alpha = 0.7) + # Set the size of the points as relative_MFI
  scale_size_continuous(
    name = "Relative MFI",
    breaks = c(0, 0.5, 1.0, 1.5),
    labels = c("0.0", "0.5", "1.0", "1.5"),
    range = c(2, 4)
  ) + 
  scale_color_manual(name = "cluster", values = c("#CDC1C5", "#96CDCD", "#DE74B0", "#FFD39B")) + # set colors
  labs(x = "Area", y = "Dispersion") + 
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white"),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank()
  ) 
## boxplot plot for relative_MFI of each state
ggplot(kmeans_result_with_cluster, aes(x = state, y = r_l)) +
  geom_boxplot(aes(color = state, fill = state), alpha = 0.7, outlier.shape = NA) + 
  scale_fill_manual(name = "cluster", values = c("#CDC1C5", "#96CDCD", "#DE74B0", "#FFD39B")) + 
  scale_color_manual(name = "cluster", values = c("#CDC1C5", "#96CDCD", "#DE74B0", "#FFD39B")) + 
  labs(x = "Cluster", y = "Relative MFI") + 
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white"),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank()
  )

# 3D plot for Area, Relative MFI, and Dispersion
# colors <- c("#CDC1C5", "#FFD39B", "#96CDCD", "#DE74B0")
# cluster_colors <- colors[kmeans_result_with_cluster$cluster]
# with(kmeans_result_with_cluster, scatter3D(x = d, y = counts, z = r_l, colvar =as.numeric(kmeans_result_with_cluster$cluster),
#                                            pch =21, cex = 1.2,col= colors,bg=cluster_colors,
#                                            xlab = "Dispersion", ylab = "Area",
#                                            zlab = "Relative MFI", 
#                                            ticktype = "detailed",bty = "b2",box = TRUE,
#                                            theta = 130, phi = 10, d=3,
#                                            colkey = F)
# )
# legend("right",title =  "classification",legend=c("no XIST", "Nebular", "Dispersed", "Confined"),pch=21,
#        cex=1.2,y.intersp=1,pt.bg = colors,bg=colors,bty="n",col=colors) 
# kmeans_result_with_cluster$state=ifelse(kmeans_result_with_cluster$cluster==1,"no XIST",ifelse(kmeans_result_with_cluster$cluster==2,"Nebular",ifelse(kmeans_result_with_cluster$cluster==3,"Dispersed",ifelse(kmeans_result_with_cluster$cluster==4,"Confined",""))))
write.table(kmeans_result_with_cluster, "XIST_kmeans_result.xls", quote=F, sep="\t", col.names=T, row.names=F)
write.table(proportion, "XIST_kmeans_proportion.xls", quote=F, sep="\t", col.names=T, row.names=F)
