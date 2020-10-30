library("clusterProfiler")
library(ggplot2)
data=read.table("../C5_mouse_epididymis_markers.tsv",sep="\t",header=T)
genes=as.character(data$gene)
cid2name=list('0'="HT1",'1'="HT2","2"="HT3","3"="HT4")
#check keyType
#library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
#mouse: org.Mm.eg.db
eg = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
name2id=list()
for(i in c(1:nrow(eg))){
  name=as.character(eg[i,1])
  id=as.character(eg[i,2])
  name2id[[name]]=id
}
gene.f=genes %in% eg$SYMBOL
genes_keep=genes[gene.f]
gene_id=sapply(genes_keep,function(x){name2id[[x]]})
gene_id=as.character(gene_id)
clusters=as.character(data$cluster)
clusters.names=sapply(clusters,function(x){cid2name[[as.character(x)]]})
my.df=data.frame(Entrez=gene_id,cluster=clusters.names[gene.f])


#for(ont_type in c("CC","MF","BP")){
ont_type="BP"
fun_type="enrichGO"
formula_res <- compareCluster(Entrez~cluster, data=my.df, fun=fun_type ,OrgDb="org.Mm.eg.db",ont=ont_type,pvalueCutoff=0.05,pAdjustMethod="BH")
d=dotplot(formula_res, x=~cluster, showCategory=6)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
pdf(paste0("Epididymis_RegionalDEGs_",fun_type,"_GO",ont_type,".pdf"),width=14,height=8)
print(d)
dev.off()
