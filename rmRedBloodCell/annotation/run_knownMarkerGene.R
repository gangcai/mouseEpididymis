library(pheatmap)
library(network)
library(ggplot2)
data=read.table("../tables/cluster_level_averageExp.tsv",header=T,sep="\t")
#marker=read.table("/home/gangcai/projects/CollaboratorLab/ChenHaoLab/mouse_epididymis/DB/Epididymis_Markers.tsv",header=T,sep="\t")
#marker=read.table("/home/gangcai/projects/CollaboratorLab/ChenHaoLab/mouse_epididymis/DB/Epididymis_Markers_Used.tsv",header=T,sep="\t")
marker=read.table("Epididymis_Markers_Used2.tsv",header=T,sep="\t")

nc=ncol(data)
mat=data[,2:nc]
mat.cname1=colnames(mat)
rownames(mat)=as.character(data[,1])
mat_marker=mat[as.character(marker[,2]),]
tag=!is.na(mat_marker[,1])
mat_marker_c=mat_marker[tag,]
group_info=as.character(apply(marker[tag,],1,function(x){paste0(x[1],"(",x[2],")")}))


pdf("marker_expression_heatmap_scaled_sampleMerged.pdf",width=6,height=7)
pheatmap(mat_marker_c,cluster_rows =F,
         scale="row",
         labels_row=group_info,
         cluster_cols = F)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
