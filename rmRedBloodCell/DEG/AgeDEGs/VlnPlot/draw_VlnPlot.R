library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(reshape2)
obj=readRDS("../../../mouse_epididymis.rds")
DefaultAssay(obj)="RNA"
clusters=obj$seurat_clusters
regions=obj$regions
cells.all=names(clusters)
ages=obj$ages

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

		cell.selected=cells.all[clusters == cid & regions == region]
		ages.selected=ages[clusters == cid & regions == region]
		metadata = data.frame("cells"=cell.selected,row.names=cell.selected,
				      "ages"=ages.selected)
		if(length(gene)>0){
			pdf(paste0("42dhigh_AgeDEGTop10_Vlnplot_",cid,"_",region,".pdf"),width=12,height=8)
			data.selected = FetchData(obj, vars=gene, cells= cell.selected, slot="data")
			obj.this <- CreateSeuratObject(t(data.selected), meta.data = metadata, 
							    project = "epididymis", min.cells = 0, min.features = 0,assay = "RNA")
		        p <- VlnPlot(obj.this,features=gene,assay =  "RNA",
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
		cell.selected=cells.all[clusters == cid & regions == region]
		ages.selected=ages[clusters == cid & regions == region]
		metadata = data.frame("cells"=cell.selected,row.names=cell.selected,
				      "ages"=ages.selected)
		if(length(gene)>0){
			pdf(paste0("56dhigh_AgeDEGTop10_Vlnplot_",cid,"_",region,".pdf"),width=12,height=8)
			data.selected = FetchData(obj, vars=gene, cells= cell.selected, slot="data")
			obj.this <- CreateSeuratObject(t(data.selected), meta.data = metadata, 
							    project = "epididymis", min.cells = 0, min.features = 0,assay = "RNA")
		        p <- VlnPlot(obj.this,features=gene,assay =  "RNA",
				      group.by="ages",pt.size=0,
		                      ncol = 5)
			print(p)
			dev.off()
		}
	}
}

