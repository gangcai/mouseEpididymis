library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
mycolor=colorRampPalette(c("tomato2","snow1","skyblue2"))(50)
obj = readRDS("../mouse_epididymis.rds")

names=list('0'="Principal cell",'1'="Myoid/Fibroblast",'2'="Clear/Narrow cell",'4'="Basal cell",
           '3'="Macrophages/Monocytes",'7'="Sperm",'6'="Endothelial cell",
           '5'="Halo/T cell")
clusters=levels(obj)

#rename the clusters
new.names <- sapply(clusters,function(x){paste0("c",x,"(",names[[as.character(x)]],")")})
obj <- RenameIdents(obj,new.names)

obj.markers <- read.table("../FindConservedMarkers/mouse_epididymis_conserved_markers_Significant.tsv",header=T,sep="\t")
#obj.markers <- obj.markers[obj.markers$pct.1 > 0.5,]
top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = logFC_min)
genes=as.character(top10$gene)
pdf("TopGenes_SingleCell_Heatmap.pdf",width=10,height=15)
p <- DoHeatmap(obj, features = genes, assay = "RNA",slot = "scale.data")+ scale_fill_gradientn(colours = rev(mycolor))
print(p)
dev.off()
