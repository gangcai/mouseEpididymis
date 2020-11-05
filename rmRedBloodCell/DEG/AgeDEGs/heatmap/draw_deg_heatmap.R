library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(pheatmap)
mycolor=colorRampPalette(c("tomato2","snow1","skyblue2"))(50)

obj=readRDS("../../../mouse_epididymis.rds")
clusters=obj$seurat_clusters
regions=obj$regions
cells.all=names(clusters)
ages=obj$ages
data.rpm=obj@assays$RNA@data #log1p(RPM)
data=read.table("../combine/AgeDEGs.tsv",header=T,sep="\t")
#Region	Cluster	Gene	Exp42d_Log1pRPM	Exp56d_Log1pRPM	Log1pRPMFC	ExpPer_42D	ExpPer_56D	Pvalue	Pvalue_Bonferroni
for(cid in c(0:7)){
	for(region in c("caput","corpus","cauda")){
		print(c(cid,region))
		filter=data$Cluster == cid & data$Region == region
		data2=data[filter,]
                data2=data2[data2$Pvalue_Bonferroni < 0.01, ]
                top10 <- data2[data2$Log1pRPMFC > 0,] %>% top_n(n = 10, wt = Log1pRPMFC)
                bottom10 <- data2[data2$Log1pRPMFC < 0, ] %>% top_n(n = -10, wt = Log1pRPMFC)
		genes_top=top10$Gene
		genes_bottom=bottom10$Gene
		genes=c(genes_top,genes_bottom)
		cells=cells.all[clusters == cid & regions == region]
		ht=5
		if(length(genes)<20){
		 ht=ht*(length(genes)/20)
		}
		pdf(paste0("AgeDEGT_Heatmap_c",cid,"_",region,".pdf"),width=10,height=ht)
		p <- DoHeatmap(obj,features=genes,cells=cells,group.by="ages",raster=F,assay="RNA",slot="scale.data")+ scale_fill_gradientn(colours = rev(mycolor))
		print(p)
		dev.off()
	}
}

