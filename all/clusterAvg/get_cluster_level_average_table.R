library(Seurat)
obj <- readRDS("../mouse_epididymis.rds")
obj.list <- SplitObject(obj, split.by = "samples")
samples=names(obj.list)

avg.df.combined=data.frame()
i=0
#get average for each sample
for(sample in samples){
	i=i+1
	obj.sample <- obj.list[[sample]]
	avg.df.sample=AverageExpression(obj.sample,assays="RNA",slot="data")
	gene=rownames(avg.df.sample$RNA)
	#write.table(cbind(gene,avg.df.sample$RNA),file=paste0(sample,"_cluster_level_averageExp.tsv"),
	#	    row.names=F,
	#	    col.names=T,sep="\t",quote=F)
	cnames=colnames(avg.df.sample$RNA)
	colnames(avg.df.sample$RNA)=sapply(cnames,function(x){paste0(sample,".c",x)})
	if(i==1){
	  avg.df.combined=avg.df.sample$RNA
	}else{
	  avg.df.combined=cbind(avg.df.combined,avg.df.sample$RNA)
	}
}

gene=rownames(avg.df.combined)
write.table(cbind(gene,avg.df.combined),file="EachSample_ClusterLevelAverage.tsv",
	    row.names=F,col.names=T,sep="\t",quote=F)

#get average for all
avg.df=AverageExpression(obj,assays="RNA",slot="data")
colnames(avg.df$RNA)=sapply(colnames(avg.df$RNA),function(x){paste0("c",x)})
gene=rownames(avg.df$RNA)
write.table(cbind(gene,avg.df$RNA),file="cluster_level_averageExp.tsv",row.names=F,
	    col.names=T,sep="\t",quote=F)
