library(pheatmap)
library(patchwork)
names=list('0'="Principal cell",'1'="Myoid/Fibroblast",'2'="Clear/Narrow cell",'4'="Basal cell",
           '3'="Macrophages/Monocytes",'7'="Sperm",'6'="Endothelial cell",
           '5'="Halo/T cell")
p.list=list()
k=0
for(cid in c(0:7)){
	data=read.table("../SegmentalDEGs_SubPopulationSpecific.tsv",header=T,sep="\t")
	filter=data[,1] == paste0("c",cid)
	data=data[filter,]
	if(nrow(data)>3){
		k=k+1
		attach(data)
		exp.df=cbind(caput_0,caput_1,caput_2,corpus_0,corpus_1,corpus_2,cauda_0,cauda_1,cauda_2)
		deg.n=nrow(exp.df)
		main.title=paste0(names[[as.character(cid)]],"(",deg.n,")")
		rownames(exp.df)=gene
		pp=pheatmap(exp.df,scale="row",
			    cluster_cols=F,cluster_rows=T,fontsize_row=5.5,
			    border_color = "white",main = main.title)
		p.list[[k]]=pp[[4]]
		detach(data)
	}
}

p=wrap_plots(p.list,ncol=4)
pdf("DEGs_heatmap.pdf",height=15,width=15)
print(p)
dev.off()
