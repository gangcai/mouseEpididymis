#epididymis Seurat scripts,
#created by: Gangcai Xie, gcxiester@gmail.com
# setting parameters

min.cell=3
min.feature=300
max.feature=4000
max.mt=50
max.ncount=20000
select.integrate.feature=3000 #not used
k.filter = 200 #not used
pc.num=30
resolution=0.1


parameters=c(min.cell,min.feature,max.feature,
             max.mt,max.ncount,select.integrate.feature,pc.num,
             resolution)
print(parameters)
suffix=paste(min.cell,min.feature,max.feature,
             max.mt,max.ncount,select.integrate.feature,pc.num,
             resolution,sep="_")
print(suffix)

library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
cluster.choose=as.character(args[1])
options(future.globals.maxSize = 30000 * 1024^2) # 30G

#load cell cycle genes
s.genes=as.character(read.table("../../../../db/cell_cycle//G1_S.mouse.txt")$V1)
g2m.genes=as.character(read.table("../../../../db/cell_cycle//G2_M.mouse.txt")$V1)
age_match=list("caput_0"="42d","caput_1"="56d","caput_2"="56d",
               "corpus_0"="42d","corpus_1"="56d","corpus_2"="56d",
               "cauda_0"="42d","cauda_1"="56d","cauda_2"="56d")

obj=readRDS("../../mouse_epididymis.rds")
#cluster.choose="0"
clusters=obj$seurat_clusters
filter=clusters == cluster.choose
counts.all=obj@assays$RNA@counts
counts.choose=counts.all[,filter]
cell.names=colnames(counts.choose)
sample.names=sapply(cell.names,function(x){strsplit(as.character(x),"[.]")[[1]][1]})
sample.names=as.character(sample.names)

regions=as.character(sapply(cell.names,function(x){ strsplit(x,"_")[[1]][1]  }))
regions=factor(regions,levels=c("caput","corpus","cauda"))

sample.names=factor(sample.names,levels=c("caput_0","caput_1","caput_2",
                                        "corpus_0","corpus_1","corpus_2",
                                        "cauda_0","cauda_1","cauda_2"))
ages=sapply(sample.names,function(x){age_match[[as.character(x)]]})
metadata = data.frame("samples"=sample.names,
                      "cells"=cell.names,
                      "ages"=ages,
                      "regions"=regions,
                      row.names=cell.names)


obj <- CreateSeuratObject(counts.choose, meta.data = metadata, project = "epididymis", min.cells = min.cell, min.features = min.feature)

obj <- PercentageFeatureSet(obj, pattern = "^mt-", col.name = "percent.mt")
p0 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "samples")
pdf(paste0("C",cluster.choose,"_epididymis_sctransform_QC_before.pdf"),height=6,width=9)
print(p0)
dev.off()


#filter
obj <- subset(obj, subset = nFeature_RNA > min.feature & nFeature_RNA < max.feature & percent.mt < max.mt & nCount_RNA<max.ncount)

obj <- SCTransform(
  obj,
  assay = 'RNA',
  new.assay.name = 'SCT',
  do.scale=T,
  do.center=T,
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA')
)
#then cacluate cell cycle scores
obj <- CellCycleScoring(
  obj,
  s.features = s.genes,
  g2m.features = g2m.genes,
  assay = 'SCT',
  set.ident = TRUE
)
#SCTranscform again by cell cycle scores and other three
obj <- SCTransform(
  obj,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score')
)


# These are now standard steps in the Seurat workflow for visualization and clustering
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:pc.num, verbose = FALSE)
obj <- RunTSNE(obj, dims = 1:pc.num, verbose = FALSE)

obj <- FindNeighbors(obj, dims = 1:pc.num, verbose = FALSE)
obj <- FindClusters(obj,resolution = resolution,  verbose = FALSE)

p <- DimPlot(obj, label = TRUE, reduction="tsne") + NoLegend()
pdf("Cluster_tSNE_Plot.pdf")
print(p)
dev.off()

p <- DimPlot(obj, label = TRUE, reduction="umap") + NoLegend()
pdf("Cluster_UMAP_Plot.pdf")
print(p)
dev.off()
#standard normalize for RNA assay
obj <- NormalizeData(obj,verbose = FALSE, normalization.method = "LogNormalize",
                                          assay="RNA",scale.factor = 1e6) #log1p(RPM)
obj <- ScaleData(obj,do.scale = TRUE , do.center = TRUE,
                                  assay="RNA")
saveRDS(obj, file =paste0("C",cluster.choose, "_mouse_epididymis.rds"))

# find markers for every cluster compared to all remaining cells, report only the positive ones
obj.markers <- FindAllMarkers(obj, assay = "RNA", slot = "data",
                                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                                 pseudocount.use = 0.5,test.use = "wilcox")
write.table(obj.markers,file=paste0("C",cluster.choose,"_mouse_epididymis_markers.tsv"),sep="\t",quote=F,row.names=F)


#QC plot
p2 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
pdf("testis_sctransform_QC_After.pdf")
print(p2)
dev.off()


#output cell ids
clusters=obj$seurat_clusters
clusters=data.frame(clusters)
write.table(clusters,file=paste0("C",cluster.choose,"_cell_cluster_ids.tsv"),sep="\t",row.names=T,col.names=F,quote=F)


#heatmap
#obj.markers <- obj.markers[obj.markers$pct.1 > 0.5,]
top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
genes=as.character(top10$gene)
pdf(paste0("C",cluster.choose,"_Top_Heatmap_EachCluster.pdf"),width=12,height=20)
p <- DoHeatmap(obj, features = genes, assay = "RNA",slot = "scale.data")
print(p)
dev.off()
