library(ggplot2)
library(patchwork)
data=read.table("../combine/AgeDEG_counts.tsv",header=T,sep="\t")
#region	cluster	clusterName	changeDirection	DEGNum
pdf("AgeDEG_Num_Barplot.pdf",width=8,height=8)
up=data[data$changeDirection == "UP",]
down=data[data$changeDirection == "DOWN",]

p1<-ggplot(up, aes(x=clusterName, y=DEGNum, fill=region)) +
  geom_bar(stat="identity",position=position_dodge())+theme_minimal() +scale_fill_brewer(palette="Dark2")+
   geom_text(aes(label=DEGNum),  vjust=-0.25, color="black",
	                 position = position_dodge(0.9), size=2)+xlab("")+ylab("# of DEGs")+
   ggtitle("56d Up Regulated")+ theme(axis.text.x = element_text(angle = 60,hjust=1,vjust=1)) 

p2<-ggplot(down, aes(x=clusterName, y=DEGNum, fill=region)) +
  geom_bar(stat="identity",position=position_dodge())+theme_minimal()+scale_fill_brewer(palette="Dark2")+
   geom_text(aes(label=DEGNum), vjust=-0.25, color="black",
	                 position = position_dodge(0.9), size=2)+xlab("")+ylab("# of DEGs")+
   ggtitle("56d Down Regulated")+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 

plots <- list(p1, p2)
wrap_plots(plots,ncol=1)
dev.off()
