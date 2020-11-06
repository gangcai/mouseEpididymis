library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
obj=readRDS("../../../../mouse_epididymis.rds")
obj$regions=factor(obj$regions,levels=c("caput","corpus","cauda"))
data=read.table("../SegmentalDEGs_SubPopulationSpecific.tsv",header=T,sep="\t")
#cluster_id      cluster_name    gene    caput_0 corpus_0        cauda_0 caput_1 corpus_1        cauda_1 caput_2 corpus_2        cauda_2 fc_0    fc_1    fc_2    padj_0  padj_1  padj_2  highestSegment  expressedPer
fc_min=apply(data[,c("fc_0","fc_1","fc_2")],1,function(x){min(x)})
padj_max=apply(data[,c("padj_0","padj_1","padj_2")],1,function(x){max(x)})
data=cbind(data,fc_min)
data=cbind(data,padj_max)
for(cid_ in c(0:7)){
	for(region in c("caput","corpus","cauda")){
		print(c(cid_,region))
		cid=paste0("c",cid_)
		filter=data$cluster_id == cid & data$highestSegment == region
		data2=data[filter,]
		#data2=data2[data2$padj_max < 0.01, ]
		top10 <- data2 %>% top_n(n = 10, wt = fc_min)
		gene=top10$gene
		if(length(gene)>10){
		 gene=gene[1:10]
		}
		ht=8
		if(length(gene)<6){
		ht=4
		}
		if(length(gene)>0){
			pdf(paste0("SegmentalDEGTop10_Vlnplot_",cid,"_",region,".pdf"),width=12,height=ht)
			p <- VlnPlot(obj,features=gene,idents=cid_,assay =  "RNA",
				group.by="regions",pt.size=0,
				ncol = 5)
			print(p)
			dev.off()
		}
	}
}


