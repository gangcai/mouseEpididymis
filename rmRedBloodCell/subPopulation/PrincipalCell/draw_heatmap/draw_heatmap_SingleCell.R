library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
mycolor=colorRampPalette(c("tomato2","snow1","skyblue2"))(50)
obj = readRDS("../C0_mouse_epididymis.rds")

names=list('0'="Prc1",'1'="Prc2","2"="Prc3","3"="Prc4","4"="Prc5","5"="Prc6","6"="Prc7","7"="Prc8")
clusters=levels(obj)

#rename the clusters
new.names <- sapply(clusters,function(x){names[[as.character(x)]]})
obj <- RenameIdents(obj,new.names)

obj.markers <- read.table("../C0_mouse_epididymis_markers.tsv",header=T,sep="\t")
top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
genes=as.character(top10$gene)
pdf("TopGenes_SingleCell_Heatmap.pdf",width=10,height=15)
p <- DoHeatmap(obj, features = genes, assay = "RNA",slot = "scale.data")+ scale_fill_gradientn(colours = rev(mycolor))
print(p)
dev.off()
