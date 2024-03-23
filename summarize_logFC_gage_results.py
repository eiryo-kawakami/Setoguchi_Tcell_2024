import math

GAGE_KO = "Setoguchi_Tcell_GAGE_KEGG_KO_ex.txt"
GAGE_Pathway = "Setoguchi_Tcell_GAGE_KEGG_Pathway_ex.txt"
GAGE_Biological_function = "Setoguchi_Tcell_GAGE_KEGG_Biological_function_ex.txt"
GAGE_Functional_catetogy = "Setoguchi_Tcell_GAGE_KEGG_Functional_category_ex.txt"

GAGE_res_list = [GAGE_KO,GAGE_Pathway,GAGE_Biological_function,GAGE_Functional_catetogy]

FuncTree_input = "Setoguchi_Tcell_GAGE_FuncTree.txt"

FDR_threshold = 0.05

with open(FuncTree_input,'w') as fo:
	for GAGE_res in GAGE_res_list:
		print GAGE_res
		with open(GAGE_res,'r') as fi:
			line = fi.readline()
			line = fi.readline()
			while line:
				itemList = line[:-1].split('\t')
				print itemList[0]
				print itemList[7]
				if float(itemList[7]) < FDR_threshold:
					#print itemList[0]
					if float(itemList[6]) > 0:
						fo.write('n-'+itemList[0].replace('\"','')+' v-'+str(-math.log10(float(itemList[7]))*50)+' o-0.6 c-#DA6272\n')
					else:
						fo.write('n-'+itemList[0].replace('\"','')+' v-'+str(-math.log10(float(itemList[7]))*50)+' o-0.6 c-#42AAC7\n')
				line = fi.readline()
