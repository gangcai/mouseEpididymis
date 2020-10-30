library(Seurat)
library(cowplot)
library(dplyr)
obj = readRDS("../C4_mouse_epididymis.rds")

names=list('0'="Bs1",'1'="Bs2","2"="Bs3","3"="Bs4")


clusters=levels(obj)

#rename the clusters
new.names <- sapply(clusters,function(x){names[[as.character(x)]]})
obj <- RenameIdents(obj,new.names)

p <- DimPlot(obj, reduction = "umap",pt.size =0.1, group.by = "samples")
pdf("Epididymis_UMAP_merged_SampleColored.pdf",width=8,height=7)
print(p)
dev.off()

p <- DimPlot(obj, label = TRUE, reduction="umap") + NoLegend()
pdf("Epididymis_UMAP_merged.pdf",width=7,height=7)
print(p)
dev.off()


obj$samples=factor(obj$samples,levels=c("caput_0","corpus_0","cauda_0","caput_1","corpus_1","cauda_1","caput_2","corpus_2","cauda_2"))
p <- DimPlot(obj, split.by="samples",reduction="umap",label = TRUE ,ncol=3) + NoLegend()
pdf("Epididymis_UMAP_each_vertical.pdf",width=12,height=12)
print(p)
dev.off()
p <- DimPlot(obj, split.by="samples",reduction="umap",label = TRUE ,ncol=3) + NoLegend()
pdf("Epididymis_UMAP_each_hondrizontal.pdf",width=12,height=12)
print(p)
dev.off()

