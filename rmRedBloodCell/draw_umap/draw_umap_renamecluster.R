library(Seurat)
library(cowplot)
library(dplyr)
obj = readRDS("../mouse_epididymis.rds")

names=list('0'="Principal cell",'1'="Myoid/Fibroblast",'2'="Clear/Narrow cell",'4'="Basal cell",
	   '3'="Macrophages/Monocytes",'7'="Sperm",'6'="Endothelial cell",
	   '5'="Halo/T cell")


clusters=levels(obj)

#rename the clusters
new.names <- sapply(clusters,function(x){paste0("c",x,"\n(",names[[as.character(x)]],")")})
obj <- RenameIdents(obj,new.names)

p <- DimPlot(obj, reduction = "umap",pt.size =0.1, group.by = "samples")
pdf("Epididymis_UMAP_merged_SampleColored.pdf",width=8,height=7)
print(p)
dev.off()

p <- DimPlot(obj, label = TRUE, reduction="umap") + NoLegend()
pdf("Epididymis_UMAP_merged.pdf",width=7,height=7)
print(p)
dev.off()

obj = readRDS("../mouse_epididymis.rds")
new.names <- sapply(clusters,function(x){paste0("c",x)})
obj <- RenameIdents(obj,new.names)

#obj$samples=factor(obj$samples,levels=c("caput","corpus","cauda"))
obj$samples=factor(obj$samples,levels=c("caput_0","corpus_0","cauda_0","caput_1","corpus_1","cauda_1","caput_2","corpus_2","cauda_2"))
p <- DimPlot(obj, split.by="samples",reduction="umap",label = TRUE ,ncol=3) + NoLegend()
pdf("Epididymis_UMAP_each_vertical.pdf",width=12,height=12)
print(p)
dev.off()
p <- DimPlot(obj, split.by="samples",reduction="umap",label = TRUE ,ncol=3) + NoLegend()
pdf("Epididymis_UMAP_each_hondrizontal.pdf",width=12,height=12)
print(p)
dev.off()

