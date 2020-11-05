import re
outfile=open("SegmentalDEGs_SubPopulationSpecific.tsv","w")

i=0
gene2info={}
for line in open("../combine/DEGs_Merged.tsv"):
	i+=1
	if i==1:
		outfile.write(line)
		continue
	items=line.rstrip().split("\t")
	cid=items[0]
	gene=items[2]
	print([gene,cid])
	if gene in gene2info:
		gene2info[gene].append(line)
	else:
		gene2info[gene]=[line]

for gene in gene2info.keys():
	lines=gene2info[gene]
	if len(lines) < 3:
		outfile.write(lines[0])
outfile.close()
