library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
obj=readRDS("../../../mouse_epididymis.rds")

data=read.table("../combine/AgeDEGs.tsv",header=T,sep="\t")
#Region	Cluster	Gene	Exp42d_Log1pRPM	Exp56d_Log1pRPM	Log1pRPMFC	ExpPer_42D	ExpPer_56D	Pvalue	Pvalue_Bonferroni
for(cid in c(0:7)){
	for(region in c("caput","corpus","cauda")){
		print(c(cid,region))
		filter=data$Cluster == cid & data$Region == region & data$Log1pRPMFC < 0
		data2=data[filter,]
		data2=data2[data2$Pvalue_Bonferroni < 0.001, ]
		top10 <- data2 %>% top_n(n = -10, wt = Log1pRPMFC)
		gene=top10$Gene
		if(length(gene)>10){
		 gene=gene[1:10]
		}
		if(length(gene)>0){
			pdf(paste0("42dhigh_AgeDEGTop10_Vlnplot_",cid,"_",region,".pdf"),width=12,height=8)
			p <- VlnPlot(obj,features=gene,idents=cid,assay =  "RNA",
				group.by="ages",pt.size=0,
				ncol = 5)
			print(p)
			dev.off()
		}
	}
}

for(cid in c(0:7)){
	for(region in c("caput","corpus","cauda")){
		print(c(cid,region))
		filter=data$Cluster == cid & data$Region == region & data$Log1pRPMFC > 0
		data2=data[filter,]
		data2=data2[data2$Pvalue_Bonferroni < 0.001, ]
		top10 <- data2 %>% top_n(n = 10, wt = Log1pRPMFC)
		gene=top10$Gene
		if(length(gene)>10){
		 gene=gene[1:10]
		}

		if(length(gene)>0){
			pdf(paste0("56dhigh_AgeDEGTop10_Vlnplot_",cid,"_",region,".pdf"),width=12,height=8)
			p <- VlnPlot(obj,features=gene,idents=cid,assay =  "RNA",
				group.by="ages",pt.size=0,
				ncol = 5)
			print(p)
			dev.off()
		}
	}
}

