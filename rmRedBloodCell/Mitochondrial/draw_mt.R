library(Seurat)
library(cowplot)
#library("RColorBrewer")
#mt.d=read.table("mt.genes.txt",header=F)
#mt.genes=mt.d$V1
#obj = readRDS("../epididymis_merged_SCT.rds")
mycolor=colorRampPalette(c("gray90","skyblue2"))(10)
obj = readRDS("../mouse_epididymis.rds")
genes=rownames(obj@assays$RNA@counts)
mt.genes=genes[grep("mt-",genes)]

names=list('0'="Principal cell",'1'="Myoid/Fibroblast",'2'="Clear/Narrow cell",'4'="Basal cell",
           '3'="Macrophages/Monocytes",'7'="Sperm",'6'="Endothelial cell",
           '5'="Halo/T cell")

clusters=levels(obj)

DefaultAssay(obj)="RNA"
#rename the clusters
new.names <- sapply(clusters,function(x){paste0("c",x,"\n(",names[[as.character(x)]],")")})
obj <- RenameIdents(obj,new.names)

p = FeaturePlot(obj, features = mt.genes,cols = c("gray90","skyblue2"),slot="data",ncol=4)
pdf("EpididymisMTGenes_FeaturePlot.pdf",height=20,width=20)
print(p)
dev.off()

p = FeaturePlot(obj, features = mt.genes,cols = mycolor,slot="data",ncol=4)
pdf("EpididymisMTGenes_FeaturePlot_color2.pdf",height=20,width=20)
print(p)
dev.off()

p = VlnPlot(obj, features = mt.genes, assay="RNA",slot="data", ncol=4)
pdf("EpididymisMTGenes_VlnPlot.pdf",height=20,width=20)
print(p)
dev.off()
