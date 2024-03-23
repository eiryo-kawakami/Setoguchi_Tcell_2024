import re

sample_info = "sample_info.txt"
sample_list = []

with open(sample_info,'r') as fi:
	line = fi.readline()
	line = fi.readline()
	while line:
		itemList = line[:-1].split('\t')
		sample_list.append(itemList[1])
		line = fi.readline()

rc = {}
gene_list = {}

for sample in sample_list:
	rc[sample] = {}
	rc_file = sample+'_rc.txt'
	with open(rc_file,'r') as fi:
		for i in range(2):
			line = fi.readline()
		while line:
			itemList = line[:-1].split('\t')
			gene = itemList[0]
			if gene not in rc[sample]:
				rc[sample][gene] = 0
			rc[sample][gene] += int(itemList[1])
			if gene not in gene_list:
				gene_list[gene] = []
			gene_list[gene].append(int(itemList[1]))
			line = fi.readline()

readcount_summary = "Setoguchi_Tcell_readcount_summary.txt"

with open(readcount_summary,'w') as fo:
	fo.write('GeneName')
	for sample in sample_list:
		fo.write('\t'+sample)
	fo.write('\n')
	for gene in gene_list:
		fo.write(gene)
		for sample in sample_list:
			fo.write('\t'+str(rc[sample][gene]))
		fo.write('\n')
