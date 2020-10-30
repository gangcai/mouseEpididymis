library(Seurat)
library(cowplot)
#library("RColorBrewer")
obj = readRDS("../mouse_epididymis.rds")
#DefaultAssay(obj)="SCT"
names=list('0'="Principal cell",'1'="Myoid/Fibroblast",'2'="Clear/Narrow cell",'4'="Basal cell",
           '3'="Macrophages/Monocytes",'7'="Sperm",'6'="Endothelial cell",
           '5'="Halo/T cell")
clusters=levels(obj)
#rename the clusters
new.names <- sapply(clusters,function(x){paste0("c",x,"\n(",names[[as.character(x)]],")")})
obj <- RenameIdents(obj,new.names)
#obj$samples=factor(obj$samples,levels=c("caput_0","corpus_0","cauda_0","caput_1","corpus_1","cauda_1","caput_2","corpus_2","cauda_2"))
obj$samples=factor(obj$samples,levels=c("caput_0","caput_1","caput_2","corpus_0","corpus_1","corpus_2","cauda_0","cauda_1","cauda_2"))
#object$samples=factor(object$samples,levels=c("caput_0","corpus_0","cauda_0","caput_1","corpus_1","cauda_1","caput_2","corpus_2","cauda_2"))
pdf("MT_per_each_population.pdf",height=6,width=15)
p <- VlnPlot(obj, features = c("percent.mt"), ncol = 1, pt.size=0, split.by = "samples")
print(p)
dev.off()
