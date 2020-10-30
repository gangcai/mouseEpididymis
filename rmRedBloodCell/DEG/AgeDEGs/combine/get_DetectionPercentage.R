library(Seurat)
library(cowplot)
library(dplyr)
library(Matrix)
obj = readRDS("../../../mouse_epididymis.rds")
counts.m=obj@assays$RNA@counts

ages=obj@meta.data$ages
regions=obj@meta.data$regions

samples=obj$samples
clusters=obj$seurat_clusters

for(cid in unique(clusters)){
	for(region in unique(regions)){
	  f.1 = clusters == cid & regions == region
	  print(c(cid,region))
	  
	  #samples.c=as.character(samples[f.1])
	  age.c=ages[f.1]
	  counts.m.c=counts.m[,f.1]
	  f1=age.c == "42d"
	  f2=age.c == "56d"
	  per=apply(counts.m.c,1,function(x){
		    x1=x[f1]
		    x2=x[f2]
		    c(sum(x1>0)/length(x1),sum(x2>0)/length(x2))
           })
	  result=cbind(rownames(counts.m.c),t(per)) 
	  colnames(result)=c("gene","42d","56d")
	  write.table(result,paste0(region,"_",cid,"_DetetionPercentage.tsv"),quote=F,sep="\t",row.names=F,col.names=T)
	}
}
