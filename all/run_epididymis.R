#epididymis Seurat scripts,
#created by: Gangcai Xie, gcxiester@gmail.com
# setting parameters
min.cell=3
min.feature=300
max.feature=4000
max.mt=50
max.ncount=20000
select.integrate.feature=3000
pc.num=50
resolution=0.05
parameters=c(min.cell,min.feature,max.feature,
	     max.mt,max.ncount,select.integrate.feature,pc.num,
	     resolution)
print(parameters)
suffix=paste(min.cell,min.feature,max.feature,
             max.mt,max.ncount,select.integrate.feature,pc.num,
             resolution,sep="_")
print(suffix)

#load libraries
library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(dplyr)
#library(doParallel)
#registerDoParallel(9)
options(future.globals.maxSize = 30000 * 1024^2) # 30G

#load cell cycle genes
s.genes=as.character(read.table("../../db/cell_cycle//G1_S.mouse.txt")$V1)
g2m.genes=as.character(read.table("../../db/cell_cycle//G2_M.mouse.txt")$V1)

age_match=list("caput_0"="42d","caput_1"="56d","caput_2"="56d",
	       "corpus_0"="42d","corpus_1"="56d","corpus_2"="56d",
	       "cauda_0"="42d","cauda_1"="56d","cauda_2"="56d")

#load mouse epididymis UMI matrix data
raw.data.merged=read.table("../../matrix/epididymis_umi_matrix.tsv",sep="\t",header=T)
all.cells=colnames(raw.data.merged)
#raw.data.merged=raw.data.merged[,all.cells=="caput_0"] #only select caput_0 for study
regions=as.character(sapply(all.cells,function(x){ strsplit(x,"_")[[1]][1]  }))
regions=factor(regions,levels=c("caput","corpus","cauda"))
all.samples=sapply(all.cells,function(x){strsplit(as.character(x),"[.]")[[1]][1]})
all.samples=as.character(all.samples)
#all.samples=factor(all.samples,levels=c("caput_0","corpus_0","cauda_0","caput_1","corpus_1","cauda_1","caput_2","corpus_2","cauda_2"))
all.samples=factor(all.samples,levels=c("caput_0","caput_1","caput_2",
					"corpus_0","corpus_1","corpus_2",
					"cauda_0","cauda_1","cauda_2"))
ages=sapply(all.samples,function(x){age_match[[as.character(x)]]})
metadata = data.frame("regions"=regions,"samples"=all.samples,
		      "cells"=all.cells,row.names=all.cells,
		      "ages"=ages)
object <- CreateSeuratObject(raw.data.merged, meta.data = metadata, project = "epididymis", min.cells = min.cell, min.features = min.feature,assay = "RNA")
#object$samples=factor(object$samples,levels=c("caput_0","corpus_0","cauda_0","caput_1","corpus_1","cauda_1","caput_2","corpus_2","cauda_2"))
object <- PercentageFeatureSet(object, pattern = "^mt-", col.name = "percent.mt")
p0 <- VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "samples")
pdf("epididymis_sctransform_QC_before.pdf",height=6,width=9)
print(p0)
dev.off()

#filter
object <- subset(object, subset = nFeature_RNA > min.feature & nFeature_RNA < max.feature & percent.mt < max.mt & nCount_RNA<max.ncount)


#SCT normalize for each sample, regress by cell cycle scores and others
object.list <- SplitObject(object, split.by = "samples")
for (i in 1:length(object.list)) {
#r <- foreach(icount(length(object.list))) %dopar% {
	#first normalized by 'percent.mt', 'nFeature_RNA', 'nCount_RNA'
        object.list[[i]] <- SCTransform(
          object.list[[i]],
          assay = 'RNA',
          new.assay.name = 'SCT',
	  do.scale=T,
	  do.center=T,
          vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA')
        )
	#then cacluate cell cycle scores
        object.list[[i]] <- CellCycleScoring(
          object.list[[i]],
          s.features = s.genes,
          g2m.features = g2m.genes,
          assay = 'SCT',
          set.ident = TRUE
        )
	#SCTranscform again by cell cycle scores and other three
        object.list[[i]] <- SCTransform(
          object.list[[i]],
          assay = 'RNA',
          new.assay.name = 'SCT',
          vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score')
        )
}
object.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = select.integrate.feature)
object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = object.features, 
				        verbose = FALSE)
object.anchors <- FindIntegrationAnchors(object.list = object.list, normalization.method = "SCT", 
					       anchor.features = object.features, verbose = FALSE)
object.integrated <- IntegrateData(anchorset = object.anchors, normalization.method = "SCT", 
				         verbose = FALSE)


# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(object.integrated) <- "integrated"

object.integrated <- RunPCA(object.integrated, npcs = pc.num, verbose = FALSE)
object.integrated <- RunUMAP(object.integrated, reduction = "pca", dims = 1:pc.num)
object.integrated <- RunTSNE(object.integrated, dims = 1:pc.num, verbose = FALSE)
object.integrated <- FindNeighbors(object.integrated, dims = 1:pc.num, verbose = FALSE)
object.integrated <- FindClusters(object.integrated,resolution = resolution,  verbose = FALSE)

#standard normalize for RNA assay
object.integrated <- NormalizeData(object.integrated,verbose = FALSE, normalization.method = "LogNormalize",
		                          assay="RNA",scale.factor = 1e6) #log1p(RPM)

object.integrated <- ScaleData(object.integrated,do.scale = TRUE , do.center = TRUE,
		                  assay="RNA")

saveRDS(object.integrated, file = "mouse_epididymis.rds")

p1 <- DimPlot(object.integrated, reduction = "umap", group.by = "samples")
p2 <- DimPlot(object.integrated, reduction = "tsne", group.by = "samples")
p3 <- DimPlot(object.integrated, label = TRUE, reduction="umap") + NoLegend()
p4 <- DimPlot(object.integrated, label = TRUE, reduction="tsne") + NoLegend()

#merged clusters with label
p5 <- DimPlot(object.integrated, split.by="samples",reduction="umap",ncol=3)
p6 <- DimPlot(object.integrated, split.by="samples",reduction="tsne",ncol=3)
#merged clusters without label
p7 <- DimPlot(object.integrated, split.by="samples",reduction="umap",label = TRUE ,ncol=3) + NoLegend()
p8 <- DimPlot(object.integrated, split.by="samples",reduction="tsne",label = TRUE, ncol=3) + NoLegend()


#QC plot
p9 <- VlnPlot(object.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "samples")

p <- plot_grid(p1, p2,p3,p4,ncol=2)
pdf(paste0("epididymis_merged_sctransform_umap_tsne_",suffix,".pdf"),height=8,width=8)
print(p)
dev.off()

pdf(paste0("epididymis_sctransform_umap_",suffix,".pdf"),height=9,width=9)
print(p5)
dev.off()

pdf(paste0("epididymis_sctransform_tsne_",suffix,".pdf"),height=9,width=9)
print(p6)
dev.off()

pdf(paste0("epididymis_sctransform_umap_labeled_",suffix,".pdf"),height=9,width=9)
print(p7)
dev.off()

pdf(paste0("epididymis_sctransform_tsne_labeled_",suffix,".pdf"),height=9,width=9)
print(p8)
dev.off()

pdf(paste0("epididymis_sctransform_QC_",suffix,".pdf"),height=6,width=9)
print(p9)
dev.off()

