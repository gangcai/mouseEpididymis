library(Seurat)
library(scsrctdb) #SingleCellSignalR based on CellTalkDB
options(future.globals.maxSize = 10000 * 1024^2) # 10G memory
obj=readRDS("../mouse_epididymis.rds")
obj.matrix=obj@assays$RNA@counts

genes.exp=rownames(obj.matrix)

clusters=obj$seurat_clusters
clusters.n=as.numeric(as.character(clusters))
#run SingleCellSignalR
cell_signal <- cell_signaling(data = obj.matrix,
                              genes = genes.exp,
                              cluster = clusters.n,
			      c.names=unique(clusters.n),
			      s.score=0.2,
			      logFC=log2(1.3),
                              gene_resive = T,
                              species = 'mus musculus')

combined.df=data.frame()
i=0
for(pair_name in names(cell_signal)){
	i=i+1
	current_signal=cell_signal[[pair_name]]
	pair_name=sub("-","-c",pair_name)
	pair_name=paste0("c",pair_name)
	this.df=cbind(rep(pair_name,nrow(current_signal)),current_signal)
	colnames(this.df)=c("cluster_pair_name","ligand","receptor","interaction_type","LRscore")
	if(i==1){
		combined.df=this.df
	}else{
		combined.df=rbind(combined.df,this.df)

	}
	write.table(current_signal,file=paste0("LR_interactions_",pair_name,"_paracrine_epididymisGenes.tsv"),sep="\t",quote=F, row.names=F)
}

write.table(combined.df,file="mouseEpididymis_ligand_receptor_communication_pairs.tsv",sep="\t",row.names=F,quote=F)

#visualize_interactions(signal = cell_signal,write.in=c(1,4))
pdf("mouse_epididymis_cell_signaling.pdf")
visualize(cell_signal)
dev.off()
for(i in names(cell_signal)){
	pdf(paste0("mouse_epididymis_cell_signaling_paracrine_show",i,".pdf"))
	visualize(cell_signal,show.in=i)
	dev.off()
}
