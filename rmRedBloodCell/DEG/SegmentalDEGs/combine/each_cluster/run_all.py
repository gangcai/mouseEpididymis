import os
for cid in range(0,8):
	print(cid)
	os.system("python combine.py {0}".format(cid))
	#os.system("Rscript draw_heatmap.R {0}".format(cid))
