library(ggplot2)
library(patchwork)
data=read.table("../DEGs_PopSpecific_Merged_count.tsv",header=T,sep="\t")
#clusterName	highestRegion	DEGNum
pdf("PopSpecificSegmentalDEG_Num_Barplot.pdf",width=6,height=4)
data$highestRegion=factor(data$highestRegion,levels=c("caput","corpus","cauda"))

p1<-ggplot(data, aes(x=clusterName, y=DEGNum, fill=highestRegion)) +
  geom_bar(stat="identity",position=position_dodge())+theme_minimal() +scale_fill_brewer(palette="Dark2")+
   geom_text(aes(label=DEGNum),  vjust=-0.25, color="black",
	                 position = position_dodge(0.9), size=2)+xlab("")+ylab("# of DEGs")+
   ggtitle("Segmental DEGs")+ theme(axis.text.x = element_text(angle = 60,hjust=1,vjust=1)) 

print(p1)
dev.off()
