library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
mycolor=colorRampPalette(c("tomato2","snow1","skyblue2"))(50)
obj = readRDS("../C4_mouse_epididymis.rds")

names=list('0'="Bs1",'1'="Bs2","2"="Bs3","3"="Bs4")
clusters=levels(obj)

#rename the clusters
new.names <- sapply(clusters,function(x){names[[as.character(x)]]})
obj <- RenameIdents(obj,new.names)

obj.markers <- read.table("../C4_mouse_epididymis_markers.tsv",header=T,sep="\t")
top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
genes=as.character(top10$gene)
pdf("TopGenes_SingleCell_Heatmap.pdf",width=10,height=15)
p <- DoHeatmap(obj, features = genes, assay = "RNA",slot = "scale.data")+ scale_fill_gradientn(colours = rev(mycolor))
print(p)
dev.off()
