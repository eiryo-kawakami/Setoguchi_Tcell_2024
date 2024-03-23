import sys, glob, os, math
import scipy.stats as stats
from rpy2.robjects import r
import rpy2.robjects as robjects

gmt_file = "/Users/eiryokawakami/Dropbox/PPI/Mouse/database/KEGG/KEGG_Pathway_gset.gmt"

gset_list = {}

with open(gmt_file,'r') as fi:
	line = fi.readline()
	while line:
		itemList = line[:-1].split('\t')
		gset_list[itemList[0]] = []
		for i in range(2,len(itemList)):
			gset_list[itemList[0]].append(itemList[i])
		line = fi.readline()

def read_genelist(genelist_file):
	gene_list = {}
	with open(genelist_file,'r') as fi:
		line = fi.readline()
		while line:
			gene_symbol = line.rstrip('\n\r').split('\t')[0]
			gene_list[gene_symbol] = 0
			line = fi.readline()
	return gene_list

def calc_fisher(positive,ind_genelist,whole_genelist):
	k = 0
	M = len(whole_genelist)
	n = 0
	N = 0
	tmp_ind_genelist = {}
	for gene in ind_genelist:
		if gene in whole_genelist:
			tmp_ind_genelist[gene] = 0
			N += 1
	for gene in positive:
		if gene in whole_genelist:
			n += 1
		if gene in tmp_ind_genelist:
			k += 1
	print "hit_in_group: " + str(k)
	print "whole_num: " + str(M)
	print "group_num: " + str(n)
	print "hit_in_whole: " + str(N)
	oddsratio, p_value = stats.fisher_exact([[k, n],[N-k,M-k-n]])
	#print p_value
	return k,N,p_value

whole_genelist_file = "Setoguchi_Tcell_logFC_gID.txt"
whole_genelist = read_genelist(whole_genelist_file)

groups = ["up","down"]

for g in groups:
	base_genes = []
	top_genes = "Setoguchi_Tcell_logFC_gID_"+g+"DEG.txt"
	with open(top_genes,'r') as fi:
		line = fi.readline()
		line = fi.readline()
		while line:
			itemList = line[:-1].split('\t')
			base_genes.append(itemList[0])
			line = fi.readline()

	p_list = []
	hit_num = []
	group_num = []
	gset_list2 = []

	for gset in gset_list:
		gset_list2.append(gset)
		k,N,p_value = calc_fisher(base_genes,gset_list[gset],whole_genelist)
		p_list.append(p_value)
		hit_num.append(k)
		group_num.append(N)

	R_p_list = robjects.FloatVector(p_list)
	r.assign('R_p_list', R_p_list)
	r.assign('list_len', len(R_p_list))
	q_list = r('p.adjust(R_p_list, method="BH", n=list_len)')

	enrichment_result = "Setoguchi_Tcell_logFC_gID_"+g+"DEG_KEGG_Pathway_enrichment.txt"

	with open(enrichment_result,'w') as fo:
		fo.write('gset\tgset_gene_num\tincluded_basis_genes\tp_value\tq_value\n')
		for j in range(len(gset_list2)):
			fo.write(gset_list2[j]+'\t')
			fo.write(str(group_num[j])+'\t')
			fo.write(str(hit_num[j])+'\t')
			fo.write(str(p_list[j])+'\t')
			fo.write(str(q_list[j])+'\n')
