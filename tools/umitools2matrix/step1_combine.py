import glob
import gzip
import re
prefix="gene"
umi_dir="./raw_from_umi_tools/"
cell2counts={}
for filename in glob.glob(umi_dir+"/*_counts.tsv.gz"):
	sample_name=filename.split("/")[-1]
	sample_name=re.sub("_counts.tsv.gz","",sample_name)
	with gzip.open(filename) as f:
		k=0
		for line in f:
			k+=1
			if k==1:
				continue
			line=line.rstrip()
			line=line.decode("utf-8")
			(transcript_id,barcode_seq,count)=line.split("\t")
			cell_id=sample_name+"."+barcode_seq
			if cell_id in cell2counts:
				cell2counts[cell_id]+=int(count)
			else:
				cell2counts[cell_id]=int(count)


cell_filter={}
min_umi=200 #minimal number of UMI
for cell_id in cell2counts.keys():
	counts=cell2counts[cell_id]
	if counts < min_umi:
		print(cell_id+" has "+str(counts)+" and total counts < "+str(min_umi))
		continue
	cell_filter[cell_id]=""


outfile=open(prefix+"_combined.tmp.tsv","w")
for filename in glob.glob(umi_dir+"/*_counts.tsv.gz"):
	sample_name=filename.split("/")[-1]
	sample_name=re.sub("_counts.tsv.gz","",sample_name)
	with gzip.open(filename) as f:
		k=0
		for line in f:
			k+=1
			if k==1:
				continue
			line=line.rstrip()
			line=line.decode("utf-8")
			(transcript_id,barcode_seq,count)=line.split("\t")
			cell_id=sample_name+"."+barcode_seq
			if cell_id in cell_filter:
				outfile.write("\t".join([transcript_id,cell_id,count])+"\n")
outfile.close()

outfile=open("cell_selected.tsv","w")
for cell_id in cell_filter.keys():
	outfile.write(cell_id+"\n")
outfile.close()
