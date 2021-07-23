import re
#caput_0_p_val   caput_0_avg_logFC       caput_0_pct.1   caput_0_pct.2   caput_0_p_val_adj       corpus_2_p_val  corpus_2_avg_logFC      corpus_2_pct.1  corpus_2_pct.2  corpus_2_p_val_adj      corpus_0_p_val  corpus_0_avg_logFC      corpus_0_pct.1  corpus_0_pct.2  corpus_0_p_val_adj      cauda_0_p_val   cauda_0_avg_logFC       cauda_0_pct.1   cauda_0_pct.2   cauda_0_p_val_adj       corpus_1_p_val  corpus_1_avg_logFC      corpus_1_pct.1  corpus_1_pct.2  corpus_1_p_val_adj      cauda_1_p_val   cauda_1_avg_logFC       cauda_1_pct.1   cauda_1_pct.2   cauda_1_p_val_adj       caput_1_p_val   caput_1_avg_logFC       caput_1_pct.1   caput_1_pct.2   caput_1_p_val_adj       cauda_2_p_val   cauda_2_avg_logFC       cauda_2_pct.1   cauda_2_pct.2   cauda_2_p_val_adj       caput_2_p_val   caput_2_avg_logFC       caput_2_pct.1   caput_2_pct.2   caput_2_p_val_adj       max_pval        minimump_p_val

def pair_hash(arr1,arr2):
	h={}
	k=-1
	for a1 in arr1:
		k+=1
		a2=arr2[k]
		h[a1]=a2
	return(h)

def get_max_padj(info_hash):
	padj_max=0
	for i in [0,1,2]:
		for region in ["caput","corpus","cauda"]:
			name=region+"_"+str(i)+"_p_val_adj"
			val=float(info_hash[name])
			if val > padj_max:
				padj_max = val
	return(str(padj_max))
def get_min_avg_logFC(info_hash):
	fc_min=0
	j=0
	for i in [0,1,2]:
		for region in ["caput","corpus","cauda"]:
			name=region+"_"+str(i)+"_avg_logFC"
			val=float(info_hash[name])
			j+=1
			if j==1:
				fc_min=val
			elif val < fc_min:
				fc_min=val
	return(str(fc_min))

#pct.1 is corresponding to the first identity
def get_mean_pct(info_hash):
	n=0
	n_sum=0
	for i in [0,1,2]:
		for region in ["caput","corpus","cauda"]:
			n+=1
			name=region+"_"+str(i)+"_pct.1"
			val=float(info_hash[name])
			n_sum+=val
	n_mean=n_sum/n
	return(str(n_mean))

def get_mean_avg_logFC(info_hash):
	n=0
	n_sum=0
	for i in [0,1,2]:
		for region in ["caput","corpus","cauda"]:
			n+=1
			name=region+"_"+str(i)+"_avg_logFC"
			val=float(info_hash[name])
			n_sum+=val
	n_mean=n_sum/n
	return(str(n_mean))


outfile=open("mouse_epididymis_conserved_markers.tsv","w")

for cid in range(0,8):
	filename="c"+str(cid)+"mouse_epididymis_markers_conserved.tsv"
	i=0
	header=""
	for line in open(filename):
		i+=1
		items=line.rstrip().split("\t")
		if i==1 and cid==0:
			header=items
			outfile.write("cluster\t"+line.rstrip()+"\tlogFC_Mean\tdetection_percentage_mean\tpvalue_adj\tlogFC_min\n")
			continue
		if i==1:
			header=items
			continue
		info_hash=pair_hash(header,items)
		logfc_mean=get_mean_avg_logFC(info_hash)
		pct_mean=get_mean_pct(info_hash)
		padj_max=get_max_padj(info_hash)
		logfc_min=get_min_avg_logFC(info_hash)
		outfile.write(str(cid)+"\t"+line.rstrip()+"\t"+logfc_mean+"\t"+pct_mean+"\t"+padj_max+"\t"+logfc_min+"\n")
outfile.close()


outfile.close()
