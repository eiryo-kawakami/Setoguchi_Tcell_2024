
log2FC = {}
pval = {}
FDR = {}
geneID = {}

pval_threshold = 0.05

logFC_file = "Setoguchi_Tcell_logFC_gID.txt"

with open(logFC_file,"r") as fi:
	line = fi.readline()
	line = fi.readline()
	while line:
		itemList = line[:-1].split('\t')
		geneID[itemList[0]] = itemList[1]
		log2FC[itemList[0]] = itemList[3]
		pval[itemList[0]] = itemList[5]
		FDR[itemList[0]] = itemList[6]
		line = fi.readline()

GO_gmt = "/Users/eiryokawakami/Library/CloudStorage/Box-Box/Per_131324_Eiryo Kawakami/PPI/Mouse/database/Uniprot/GO_gset_biological_process.gmt"

GO = {}
GO_name = {}

with open(GO_gmt,"r") as fi:
	line = fi.readline()
	while line:
		itemList = line[:-1].split('\t')
		GO[itemList[1]] = []
		GO_name[itemList[1]] = itemList[0]
		for i in range(2,len(itemList)):
			GO[itemList[1]].append(itemList[i])
		line = fi.readline()

GO_Enrichment = "Setoguchi_Tcell_GAGE_GO_biological_process.txt"

up_GO = []
down_GO = []

with open(GO_Enrichment,'r') as fi:
	line =fi.readline()
	line =fi.readline()
	while line:
		itemList = line[:-1].split('\t')
		if float(itemList[3]) < pval_threshold:
			if float(itemList[2]) > 0:
				up_GO.append(itemList[0].split(' ')[0].replace('\"',''))
			else:
				down_GO.append(itemList[0].split(' ')[0].replace('\"',''))
		line = fi.readline()

annt = logFC_file.replace('.txt',"_annotated.txt")
with open(logFC_file,"r") as fi:
	with open(annt,"w") as fo:
		fo.write("GeneName\tGeneID\tlog2FC\tpval\tFDR")
		for p in up_GO:
			fo.write('\t'+GO_name[p])
		for p in down_GO:
			fo.write('\t'+GO_name[p])
		fo.write('\n')
		line = fi.readline()
		line = fi.readline()
		while line:
			gene = line[:-1].split('\t')[0]
			fo.write(gene+'\t'+geneID[gene]+'\t'+log2FC[gene]+'\t'+pval[gene]+'\t'+FDR[gene])
			for p in up_GO:
				if geneID[gene] in GO[p]:
					fo.write('\t1')
				else:
					fo.write('\t0')
			for p in down_GO:
				if geneID[gene] in GO[p]:
					fo.write('\t1')
				else:
					fo.write('\t0')
			fo.write('\n')
			line = fi.readline()




