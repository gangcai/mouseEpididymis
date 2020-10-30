import re,sys
expressed_per_cutoff=0.25 #cutoff for batch minimal and regional maximal percentage of expressed genes
fc_cutoff=2 # log((1+RPM1)/(1+RPM2)) change cutoff
pvalue_cutoff=0.05 #batch maximal p value (adjusted for multiple comparison) should be less than this
cid=sys.argv[1]
sample2info={}
all_genes={}
for sample in ["0","1","2"]:
	filename="../../Sample"+sample+"_cluster"+cid+"_KruskalTest.tsv"
	i=0
	sample2info[sample]={}
	for line in open(filename):
		i+=1
		if i==1:
			continue
		items=line.rstrip().split()
		gene=items[0]
		all_genes[gene]=""
		fc=items[1]
		padj=items[-1]
		pvalue=items[2]
		caput=items[3]
		cauda=items[4]
		corpus=items[5]
		sample2info[sample][gene]=[caput,corpus,cauda,fc,pvalue,padj]


gene2per={} #gene id to the regional maximal percentage of cells expressed at each batch
gene2per_detail={}
gene2per_detail2={}
for sample in ["0","1","2"]:
	if not sample in gene2per_detail:
		gene2per_detail[sample]={}
	for region in ["caput","corpus","cauda"]:
		filename="../../GetDetectionProportion/"+region+"_"+sample+"_"+cid+"_DetetionPercentage.tsv"
		i=0
		for line in open(filename):
			i+=1
			if i==1:
				continue
			(gene,per)=line.rstrip().split("\t")
			if not gene in gene2per_detail[sample]:
				gene2per_detail[sample][gene]=[float(per)]
			else:
				gene2per_detail[sample][gene].append(float(per))
print(len(gene2per_detail["0"].keys()))
for sample in gene2per_detail.keys():
	for gene in gene2per_detail[sample].keys():
		per_max=max(gene2per_detail[sample][gene])
		if not gene in gene2per_detail2:
			gene2per_detail2[gene]=[per_max]
		else:
			gene2per_detail2[gene].append(per_max)

for gene in gene2per_detail2.keys():
	per=min(gene2per_detail2[gene]) #minimal detection percentage among the regional maximal percentage of all batches
	gene2per[gene]=per



def get_max_region_id(caput,corpus,cauda):
	caput=float(caput)
	corpus=float(corpus)
	cauda=float(cauda)
	if caput > corpus and caput > cauda:
		return("caput")

	if corpus > caput and corpus > cauda:
		return("corpus")

	if cauda > caput and cauda > corpus:
		return("cauda")

outfile=open("cluster"+cid+"_tested_info_combined.tsv","w")
outfile.write("gene\tcaput_0\tcorpus_0\tcauda_0\tcaput_1\tcorpus_1\tcauda_1\tcaput_2\tcorpus_2\tcauda_2\tfc_0\tfc_1\tfc_2\tpvalue_0\tpvalue_1\tpvalue_2\tpadj_0\tpadj_1\tpadj_2\thighestSegment_0\thighestSegment_1\thighestSegment_2\texpressedPer\n")

outfile2=open("cluster"+cid+"_testedSignificant_info_combined.tsv","w")
outfile2.write("gene\tcaput_0\tcorpus_0\tcauda_0\tcaput_1\tcorpus_1\tcauda_1\tcaput_2\tcorpus_2\tcauda_2\tfc_0\tfc_1\tfc_2\tpvalue_0\tpvalue_1\tpvalue_2\tpadj_0\tpadj_1\tpadj_2\thighestSegment\texpressedPer\n")
def check_sig_fc(fcs):
	if "NA" in fcs:
		return(False)

	min_fc=min([float(k) for k in fcs])
	if min_fc > fc_cutoff:
		return(True)
	else:
		return(False)

def check_consistency(max_regions):
	if "NA" in max_regions:
		return(False)
	unique={}
	for mr in max_regions:
		unique[mr]=""
	if len(unique.keys())==1:
		return(True)
	else:
		return(False)

def check_pvalues(pvalues):
	if "NA" in pvalues:
		return(False)
	max_pvalue=max([float(k) for k in pvalues])
	if max_pvalue < pvalue_cutoff:
		return(True)
	else:
		return(False)
	

for gene in all_genes.keys():
	exps=[]
	fcs=[]
	padjs=[]
	pvalues=[]
	max_regions=[]
	for sample in ["0","1","2"]:
		if gene in sample2info[sample]:
			(caput,corpus,cauda,fc,pvalue,padj)=sample2info[sample][gene]
			exps.append(caput)
			exps.append(corpus)
			exps.append(cauda)
			fcs.append(fc)
			pvalues.append(pvalue)
			padjs.append(padj)
			max_region=get_max_region_id(caput,corpus,cauda)
			max_regions.append(max_region)
		else:
			(caput,corpus,cauda,fc,pvalue,padj)=["NA","NA","NA","NA","NA","NA"]
			exps.append(caput)
			exps.append(corpus)
			exps.append(cauda)
			fcs.append(fc)
			padjs.append(padj)
			pvalues.append(pvalue)
			max_regions.append("NA")
	expressed_per=gene2per[gene]
	if check_sig_fc(fcs) and check_consistency(max_regions) and check_pvalues(padjs) and expressed_per > expressed_per_cutoff:
		outfile2.write(gene+"\t"+"\t".join(exps)+"\t"+"\t".join(fcs)+"\t"+"\t".join(pvalues)+"\t"+"\t".join(padjs)+"\t"+max_regions[0]+"\t"+str(expressed_per)+"\n")
	outfile.write(gene+"\t"+"\t".join(exps)+"\t"+"\t".join(fcs)+"\t"+"\t".join(pvalues)+"\t"+"\t".join(padjs)+"\t"+"\t".join(max_regions)+"\t"+str(expressed_per)+"\n")
outfile.close()
outfile2.close()
