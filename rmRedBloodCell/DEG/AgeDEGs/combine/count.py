import re
names={'0':"Principal cell",'1':"Myoid/Fibroblast",'2':"Clear/Narrow cell",'4':"Basal cell",
		           '3':"Macrophages/Monocytes",'7':"Sperm",'6':"Endothelial cell",
			              '5':"Halo/T cell"}
outfile=open("AgeDEG_counts.tsv","w")
i=0
count={}
for line in open("AgeDEGs.tsv"):
	i+=1
	if i==1:
		continue
	items=line.rstrip().split()
	(region,cluster)=items[0:2]
	logfc=float(items[-5])
	if logfc > 0:
		change="UP" #higher in 56d
	else:
		change="DOWN" #lower in 56d
	cluster_name=names[cluster]
	id_=region+"\t"+cluster+"\t"+cluster_name+"\t"+change
	if id_ in count:
		count[id_]+=1
	else:
		count[id_]=1
outfile.write("region\tcluster\tclusterName\tchangeDirection\tDEGNum\n")
for id_ in count.keys():
	num=str(count[id_])
	outfile.write(id_+"\t"+num+"\n")

outfile.close()
