library(Seurat)
library(cowplot)
library(dplyr)
library(Matrix)
#options <- commandArgs(trailingOnly = TRUE)
obj = readRDS("../../mouse_epididymis.rds")
counts.m=obj@assays$RNA@data #log1p(RPM)

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
	  print(dim(counts.m.c))
	  
	  
	  #kruskal test
	  kt.test <- apply(counts.m.c,1,function(x){
	    dt = data.frame(exp=as.numeric(x),age=age.c)
	    kt = kruskal.test(exp~age,data=dt)
	    pvalue = kt$p.value
	    exp.avg = aggregate(dt[,"exp"], list(dt$age), mean)$x
	    exp.max=max(exp.avg)
	    exp.min=min(exp.avg)
	    fc=exp.max - exp.min # log((1+RPM1)/(1+RPM2)) = log1pRPM1 - log1pRPM2
	    #fc=(exp.max+0.1)/(exp.min+0.1)
	    c(fc,pvalue,exp.avg)
	  })
	  dt = data.frame(exp=as.numeric(counts.m.c[1,]),age=age.c)
	  exp.avg.name=as.character(aggregate(dt[,"exp"], list(dt$age), mean)$Group.1)
	  kt.test=t(kt.test)
	  
	  p.filter=!is.na(kt.test[,2])
	  
	  kt.test.2=kt.test[p.filter,]
	  p.adj=p.adjust(as.numeric(kt.test.2[,2]),method="bonferroni")
	  
	  genes=rownames(kt.test.2)
	  
	  results=cbind(genes,kt.test.2,p.adj)
	  
	  pa.filter=is.na(p.adj)
	  
	  
	  colnames(results)=c("gene","foldchange","pvalue",exp.avg.name,"pvalue.bonferroni")
	  write.table(results,file=paste0(region,"_cluster",cid,"_KruskalTest.tsv"),sep="\t",
		      col.names = T,row.names = F,
		      quote=F)
	}
}
