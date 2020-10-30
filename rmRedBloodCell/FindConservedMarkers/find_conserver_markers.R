library(Seurat)
#options(future.globals.maxSize = 600 * 1024^2)
args = commandArgs(trailingOnly=TRUE)
cid=as.character(args[1])
object.integrated=readRDS("../mouse_epididymis.rds")
DefaultAssay(object.integrated) <- "RNA"
#get marker genes for each cluster
object.markers <- FindConservedMarkers(object.integrated,ident.1=cid,grouping.var="samples")
gene=rownames(object.markers)
marker.genes=cbind(gene,object.markers)
filename=paste0("c",cid,"mouse_epididymis_markers_conserved.tsv")
write.table(marker.genes,file=filename,sep="\t",quote=F,row.names=F)
