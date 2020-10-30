import glob,re
k=0
log2fc_cutoff=1 # log2fc= log2((1+rpm1)/(1+rpm2))
cutoff=0.01
percentage_cutoff=0.5 # maximal proportion of cell detected in either age group should be larger than this
def gene2per(fname):
	j=0
	h={}
	for line in open(fname):
		j+=1
		if j==1:
			continue
		(gene,per1,per2)=line.rstrip().split()
		h[gene]=[per1,per2]
	return(h)

outfile=open("AgeDEGs.tsv","w")
for filename in glob.glob("../*_KruskalTest.tsv"):
	fname=filename.split("/")[-1]
	(region,cluster,other)=fname.split("_")

	f2=region+"_"+re.sub("cluster","",cluster)+"_DetetionPercentage.tsv"
	g2per=gene2per(f2)
	k+=1
	i=0
	#gene    foldchange      pvalue  42d     56d     pvalue.bonferroni
	for line in open(filename):
		i+=1
		if i==1:
			if k==1:
				#outfile.write("region\tcluster\t"+line.rstrip()+"\t42d_Per\t56d_per\n")
				outfile.write("Region\tCluster\tGene\tExp42d_Log1pRPM\tExp56d_Log1pRPM\tLog1pRPMFC\tExpPer_42D\tExpPer_56D\tPvalue\tPvalue_Bonferroni\n")
			continue
		#items=line.rstrip().split("\t")
		(gene,foldchange,pvalue,exp_42d,exp_56d,p_adj)=line.rstrip().split("\t")
		#log1p_fc=float(items[-2])-float(items[-1])
		log1p_fc=float(exp_56d)-float(exp_42d)
		if abs(log1p_fc) < log2fc_cutoff:
			continue
		(per1,per2)=g2per[gene]
		per_max=max([float(per1),float(per2)])
		#Region  Cluster Gene    Exp42d_Log1pRPM Exp56d_Log1pRPM Log1pRPMFC      ExpPer_42D      ExpPer_56D
		results=[region,re.sub("cluster","",cluster),gene,exp_42d,exp_56d,str(log1p_fc),per1,per2,pvalue,p_adj]
		if float(p_adj) < cutoff and per_max > percentage_cutoff:
			cluster=re.sub("cluster","",cluster)
			#outfile.write(region+"\t"+cluster+"\t"+line.rstrip()+"\t"+per1+"\t"+per2+"\n")
			outfile.write("\t".join(results)+"\n")
			continue
outfile.close()
