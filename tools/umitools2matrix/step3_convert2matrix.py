sample_names_d={}
for line in open("cell_selected.tsv"):
	cell_name=line.rstrip()
	sample_names_d[cell_name]=""

samples=sample_names_d.keys()
samples=sorted(samples)

outfile=open("Gene_Exp_Matrix.tsv","w")
outfile.write("\t".join(samples)+"\n") #no first column header
pre_name=""
i=0
pre_hash={}
exp_cutoff=0 # larger than 0 is defined as expressed
min_num_of_cells=3 #at least this number of cell has expression
#only output the genes expressed in at least 3 cell
for line in open("gene_combined_sorted.tmp.tsv"):
	line=line.rstrip()
	i+=1
	(gene_id,cell_name,count)=line.split("\t")
	if i==1:
		pre_name=gene_id
		pre_hash[cell_name]=count
		continue
	if pre_name == gene_id:
		pre_hash[cell_name]=count
		continue

	f_cell_num=0 #filtered cell number
	for c_n in pre_hash.values():
		c_n=float(c_n)
		if c_n > exp_cutoff:
			f_cell_num+=1 # at least count 1

	if f_cell_num >= min_num_of_cells: 
		count_a=[pre_name]
		for sname in samples:
			if sname in pre_hash:
				this_count=pre_hash[sname]
			else:
				this_count="0"
			count_a.append(this_count)
		outfile.write("\t".join(count_a)+"\n")
	#reset
	pre_name=gene_id
	pre_hash={}
	pre_hash[cell_name]=count

#the last one
f_cell_num=0 #filtered cell number
for c_n in pre_hash.values():
	c_n=float(c_n)
	if c_n > exp_cutoff:
		f_cell_num+=1 # at least count 1

if f_cell_num >= min_num_of_cells: 
	count_a=[pre_name]
	for sname in samples:
		if sname in pre_hash:
			this_count=pre_hash[sname]
		else:
			this_count="0"
		count_a.append(this_count)
	outfile.write("\t".join(count_a)+"\n")
outfile.close()
