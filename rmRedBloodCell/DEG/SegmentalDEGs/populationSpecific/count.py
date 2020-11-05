import re
names={'0':"Principal cell",'1':"Myoid/Fibroblast",'2':"Clear/Narrow cell",'4':"Basal cell",
		           '3':"Macrophages/Monocytes",'7':"Sperm",'6':"Endothelial cell",
			              '5':"Halo/T cell"}

count={}
for cid in names.values():
	for seg in ["caput","corpus","cauda"]:
		id_=cid+"\t"+seg
		count[id_]=0


outfile=open("DEGs_PopSpecific_Merged_count.tsv","w")
i=0
for line in open("SegmentalDEGs_SubPopulationSpecific.tsv"):
	i+=1
	if i==1:
		continue
	items=line.rstrip().split("\t")
	highestSegment=items[-2]
	cluster_name=items[1]
	id_=cluster_name+"\t"+highestSegment
	count[id_]+=1
outfile.write("clusterName\thighestRegion\tDEGNum\n")
for id_ in count.keys():
	num=str(count[id_])
	outfile.write(id_+"\t"+num+"\n")

outfile.close()
