# Setoguchi_Tcell_2024
Code and intermediate data used in Setoguchi et al. Nat Comm 2024.

**execSTAR.sh**
[Reference] UCSC mm10
[Input] [samplename]_H1_L001_R1.fastq.gz
[Output] [samplename].Aligned.sortedByCoord.out.bam

**Rsubread.R**
[Reference] UCSC mm10 genes.gtf
[Input] [samplename].Aligned.sortedByCoord.out.bam
[Output] [samplename]_rc.txt, [samplename]_rpkm.txt

**summarize_rc.py**
[Input] sample_info.txt, [samplename]_rc.txt
[Output] Setoguchi_Tcell_readcount_summary.txt

**execDESeq2.txt**
[Input] sample_info.txt, Setoguchi_Tcell_readcount_summary.txt
[Output] Setoguchi_Tcell_logFC.txt, Setoguchi_Tcell_rlog_count.txt, Setoguchi_Tcell_rlog_count_center_cutoff.txt

**PCA_all.R**
[Input] sample_info.txt Setoguchi_Tcell_rlog_count_center_cutoff.txt
[Output] Setoguchi_Tcell_PCA.pdf, Setoguchi_Tcell_tSNE.pdf

**sym2gID.R**
[Input] Setoguchi_Tcell_logFC.txt, Setoguchi_Tcell_rlog_count_center_cutoff.txt
[Output] Setoguchi_Tcell_logFC_gID.txt, Setoguchi_Tcell_rlog_count_center_cutoff_gID.txt

**volcanoPlot.R**
[Input] Setoguchi_Tcell_logFC_gID.txt
[Output] Setoguchi_Tcell_logFC_gID_upDEG.txt, Setoguchi_Tcell_logFC_gID_downDEG.txt, Setoguchi_Tcell_gID_volcanoPlot.pdf

**KEGGenrichment.py**
[Reference] KEGG_Pathway_gset.gmt
[Input] Setoguchi_Tcell_logFC_gID.txt, Setoguchi_Tcell_logFC_gID_upDEG.txt, Setoguchi_Tcell_logFC_gID_downDEG.txt
[Output] Setoguchi_Tcell_logFC_gID_upDEG_KEGG_Pathway_enrichment.txt, Setoguchi_Tcell_logFC_gID_downDEG_KEGG_Pathway_enrichment.txt

**GO_GAGE_original_gmt2.R**
[Reference] GO_gset_biological_process.gmt, GO_gset_cellular_component.gmt, GO_gset_molecular_function.gmt
[Input] Setoguchi_Tcell_logFC_gID.txt
[Output] Setoguchi_Tcell_GAGE_GO_biological_process.txt, Setoguchi_Tcell_GAGE_GO_cellular_component.txt, Setoguchi_Tcell_GAGE_GO_molecular_function.txt

**KEGG_GAGE.R**
[Reference] KEGG_KO_gset.gmt, KEGG_Pathway_gset.gmt, KEGG_Biological_function_gset.gmt, KEGG_Functional_category_gset.gmt
[Input] Setoguchi_Tcell_logFC_gID.txt
[Output] Setoguchi_Tcell_GAGE_KEGG_KO.txt, Setoguchi_Tcell_GAGE_KEGG_Pathway.txt, Setoguchi_Tcell_GAGE_KEGG_Biological_function.txt

**KEGG_lm.R**
[Input] sample_info.txt, Setoguchi_Tcell_GAGE_KEGG_Biological_function.txt, Setoguchi_Tcell_GAGE_KEGG_Functional_category.txt, Setoguchi_Tcell_GAGE_KEGG_Pathway.txt, Setoguchi_Tcell_GAGE_KEGG_KO.txt
[Output] Setoguchi_Tcell_GAGE_KEGG_Biological_function_ex.txt, Setoguchi_Tcell_GAGE_KEGG_Functional_category_ex.txt, Setoguchi_Tcell_GAGE_KEGG_Pathway_ex.txt, Setoguchi_Tcell_GAGE_KEGG_KO_ex.txt

**summarize_gage_results.py**
[Input] Setoguchi_Tcell_GAGE_KEGG_KO_ex.txt, Setoguchi_Tcell_GAGE_KEGG_Pathway_ex.txt, Setoguchi_Tcell_GAGE_KEGG_Biological_function_ex.txt, Setoguchi_Tcell_GAGE_KEGG_Functional_category_ex.txt
[Output] Setoguchi_Tcell_GAGE_FuncTree.txt

**selectGO.py**
[Input] Setoguchi_Tcell_GAGE_GO_biological_process.txt
[Output] Setoguchi_Tcello_GAGE_GO_biological_process_down_GOvis.txt

**wPGSA.py**
[Usage] python wPGSA.py --network-file merged_mouse_v7_160212.network --logfc-file Setoguchi_Tcell_logFC.txt
[Input] Setoguchi_Tcell_logFC.txt
[Output] Setoguchi_Tcell_logFC_TF_wPGSA_p_value.txt, Setoguchi_Tcell_logFC_TF_wPGSA_q_value.txt, Setoguchi_Tcell_logFC_TF_wPGSA_t_score.txt

**annotate_genes.py**
[Reference] GO_gset_biological_process.gmt
[Input] Setoguchi_Tcell_logFC_gID.txt, Setoguchi_Tcell_GAGE_GO_biological_process.txt
[Output] Setoguchi_Tcell_GAGE_GO_biological_process_annotated.txt

**pathview.R**
[Input] Setoguchi_Tcell_logFC_gID.txt, Setoguchi_Tcell_GAGE_KEGG_Pathway_ex.txt
[Output] mmuxxxxx.png, mmuxxxxx..multi.png

