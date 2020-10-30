library(Seurat)
library(cowplot)
library(dplyr)
library(Matrix)
obj = readRDS("../../../mouse_epididymis.rds")
counts.m=obj@assays$RNA@counts #UMI
samples=obj$samples
clusters=obj$seurat_clusters

for(cid in unique(clusters)){
	for(sample in unique(samples)){
	  f.1 = clusters == cid & samples == sample
	  print(c(cid,sample))
	  counts.m.c=counts.m[,f.1]
	  per=apply(counts.m.c,1,function(x){
		    sum(x>0)/length(x)
           })
	  result=cbind(rownames(counts.m.c),as.numeric(per)) 
	  colnames(result)=c("gene","Percentage")
	  write.table(result,paste0(sample,"_",cid,"_DetetionPercentage.tsv"),quote=F,sep="\t",row.names=F,col.names=T)
	}
}
