library(Seurat)
library(cowplot)
#library("RColorBrewer")
library(dplyr)
library(ggplot2)
#mycolor=colorRampPalette(c("tomato2","snow1","skyblue2"))(10)
mycolor=colorRampPalette(c("skyblue2","snow1","tomato2"))(10)
obj = readRDS("../C0_mouse_epididymis.rds")

names=list('0'="Prc1",'1'="Prc2","2"="Prc3","3"="Prc4","4"="Prc5","5"="Prc6","6"="Prc7","7"="Prc8")


clusters=levels(obj)

#rename the clusters
new.names <- sapply(clusters,function(x){names[[as.character(x)]]})
obj <- RenameIdents(obj,new.names)

DefaultAssay(obj) <- "RNA"

p = FeaturePlot(obj, features = c("Actb"), cols = mycolor,slot="data",ncol=1)
pdf("Actb_PrincipalSubpop_FeaturePlot.pdf",height=5,width=5)
print(p)
dev.off()

p = VlnPlot(obj,assay = "RNA", features = c("Actb"), slot="data", ncol=1,pt.size=0)
pdf("Actb_PrincipalSubpop_VlnPlot.pdf",height=8,width=6)
print(p)
dev.off()



obj.markers <- read.table("../C0_mouse_epididymis_markers.tsv",header=T,sep="\t")
top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
gene=top10[top10$cluster == "6",]$gene
p = FeaturePlot(obj, features = gene, cols = mycolor,slot="data",ncol=5)
pdf("Prc7_PrincipalSubpop_FeaturePlot.pdf",height=6,width=15)
print(p)
dev.off()

p = VlnPlot(obj,assay = "RNA", features = gene, slot="data", ncol=5,pt.size=0)
pdf("Prc7_PrincipalSubpop_VlnPlot.pdf",height=8,width=15)
print(p)
dev.off()


p = FeaturePlot(obj, features = gene, cols = mycolor,slot="data",ncol=2)
pdf("Prc7_PrincipalSubpop_FeaturePlot_v2.pdf",height=15,width=6)
print(p)
dev.off()

p = VlnPlot(obj,assay = "RNA", features = gene, slot="data", ncol=2,pt.size=0)
pdf("Prc7_PrincipalSubpop_VlnPlot_v2.pdf",height=15,width=6)
print(p)
dev.off()


