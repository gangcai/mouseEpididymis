data=read.table("mouse_epididymis_conserved_markers.tsv",header=T,sep="\t")
filter=data$pvalue_adj  < 0.05
write.table(data[filter,],file="mouse_epididymis_conserved_markers_Significant.tsv",quote=F,row.names=F,col.names=T,sep="\t")
