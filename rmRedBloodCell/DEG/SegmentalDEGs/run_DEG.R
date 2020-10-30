library(Seurat)
library(cowplot)
library(dplyr)
library(Matrix)
options <- commandArgs(trailingOnly = TRUE)
batch=as.character(options[1])
obj = readRDS("../../mouse_epididymis.rds")
counts.m.all=obj@assays$RNA@data #log1p(RPM)
cell_names=colnames(counts.m.all)
cell_names_choose=cell_names[grep(batch,cell_names)]
counts.m=counts.m.all[,cell_names_choose]

samples.all=obj$samples
clusters.all=obj$seurat_clusters
samples.grep=grep(batch,samples.all)
clusters=clusters.all[samples.grep] #choose given batch
samples=samples.all[samples.grep]

for(cid in unique(clusters)){
  
  f.1 = clusters == cid
  print(cid)
  
  samples.c=as.character(samples[f.1])
  
  counts.m.c=counts.m[,f.1]
  print(dim(counts.m.c))
  
  
  #kruskal test
  kt.test <- apply(counts.m.c,1,function(x){
    dt = data.frame(exp=as.numeric(x),region=samples.c)
    kt = kruskal.test(exp~region,data=dt)
    pvalue = kt$p.value
    exp.avg = aggregate(dt[,"exp"], list(dt$region), mean)$x
    exp.max=max(exp.avg)
    exp.min=min(exp.avg)
    #fc=(exp.max+0.1)/(exp.min+0.1)
    fc=exp.max-exp.min  #log((1+RPM1)/(1+RPM2)) = log1pRPM1 - log1pRPM2
    c(fc,pvalue,exp.avg)
  })
  dt = data.frame(exp=as.numeric(counts.m.c[1,]),region=samples.c)
  exp.avg.name=as.character(aggregate(dt[,"exp"], list(dt$region), mean)$Group.1)
  kt.test=t(kt.test)
  
  p.filter=!is.na(kt.test[,2])
  
  kt.test.2=kt.test[p.filter,]
  p.adj=p.adjust(as.numeric(kt.test.2[,2]),method="bonferroni")
  
  genes=rownames(kt.test.2)
  
  results=cbind(genes,kt.test.2,p.adj)
  
  pa.filter=is.na(p.adj)
  
  
  colnames(results)=c("gene","foldchange","pvalue",exp.avg.name,"pvalue.bonferroni")
  write.table(results,file=paste0("Sample",batch,"_cluster",cid,"_KruskalTest.tsv"),sep="\t",
              col.names = T,row.names = F,
              quote=F)
}

