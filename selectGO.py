GO_summary = "Setoguchi_Tcello_GAGE_GO_biological_process_down_GOvis.txt"

pval_threshold = 0.05

up_color = "#BC2638"
down_color = "#325CA0"

GOenrichment = "Setoguchi_Tcell_GAGE_GO_biological_process.txt"

up_GO = []
down_GO = []

with open(GOenrichment,'r') as fi:
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

with open(GO_summary,'w') as fo:
	for GO in up_GO:
		fo.write(GO+' '+up_color+'\n')
	for GO in down_GO:
		fo.write(GO+' '+down_color+'\n')


