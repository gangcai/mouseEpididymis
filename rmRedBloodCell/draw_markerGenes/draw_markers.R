library(Seurat)
library(cowplot)
#library("RColorBrewer")
library(dplyr)
library(ggplot2)
#mycolor=colorRampPalette(c("tomato2","snow1","skyblue2"))(10)
mycolor=colorRampPalette(c("skyblue2","snow1","tomato2"))(10)
obj = readRDS("../mouse_epididymis.rds")
names=list('0'="Principal cell",'1'="Myoid/Fibroblast",'2'="Clear/Narrow cell",'4'="Basal cell",
           '3'="Macrophages/Monocytes",'7'="Sperm",'6'="Endothelial cell",
           '5'="Halo/T cell")
clusters=levels(obj)
DefaultAssay(obj) <- "RNA"
#rename the clusters
new.names <- sapply(clusters,function(x){names[[as.character(x)]]})
obj <- RenameIdents(obj,new.names)

p = FeaturePlot(obj, features = c("Krt14","C4bp","Adam7"), cols = mycolor,slot="data",ncol=1)
pdf("Aqp9_EpididymisCellMarker_FeaturePlot.pdf",height=5,width=15)
print(p)
dev.off()

p = VlnPlot(obj,assay = "RNA", features = c("Krt14","C4bp","Adam7"), slot="data", ncol=3,pt.size=0.5)
pdf("EpididymisCellMarker_VlnPlot.pdf",height=5,width=15)
print(p)
dev.off()

p = VlnPlot(obj,assay = "RNA", features = c("Krt14","C4bp","Adam7"), slot="data", ncol=3,pt.size=0)
pdf("EpididymisCellMarker_VlnPlot_ptsize0.pdf",height=5,width=15)
print(p)
dev.off()

p = FeaturePlot(obj, features = c("Aqp9"), cols = mycolor,slot="data",ncol=1)
pdf("Aqp9_EpididymisCellMarker_FeaturePlot.pdf",height=5,width=5)
print(p)
dev.off()

p = VlnPlot(obj,assay = "RNA", features = c("Aqp9"), slot="data", ncol=1, pt.size=0.5)
pdf("Aqp9_EpididymisCellMarker_VlnPlot.pdf",height=5,width=5)
print(p)
dev.off()

p = FeaturePlot(obj, features = c("Krt5"), cols = mycolor,slot="data",ncol=1)
pdf("Krt5_EpididymisCellMarker_FeaturePlot.pdf",height=5,width=5)
print(p)
dev.off()

p = VlnPlot(obj,assay = "RNA", features = c("Krt5"), slot="data", ncol=1, pt.size=0.5)
pdf("Krt5_EpididymisCellMarker_VlnPlot.pdf",height=5,width=5)
print(p)
dev.off()

p = FeaturePlot(obj, features = c("Cd68"), cols = mycolor,slot="data",ncol=1)
pdf("Cd68_EpididymisCellMarker_FeaturePlot.pdf",height=5,width=5)
print(p)
dev.off()

p = VlnPlot(obj,assay = "RNA", features = c("Cd68"), slot="data", ncol=1, pt.size=0.5)
pdf("Cd68_EpididymisCellMarker_VlnPlot.pdf",height=5,width=5)
print(p)
dev.off()


p = FeaturePlot(obj, features = c("Atp6v1b1"), cols = mycolor,slot="data",ncol=1)
pdf("Atp6v1b1_EpididymisCellMarker_FeaturePlot.pdf",height=5,width=5)
print(p)
dev.off()

p = VlnPlot(obj,assay = "RNA", features = c("Atp6v1b1"), slot="data", ncol=1, pt.size=0.5)
pdf("Atp6v1b1_EpididymisCellMarker_VlnPlot.pdf",height=5,width=5)
print(p)
dev.off()
