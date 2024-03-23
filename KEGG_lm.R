sample_info <- read.table("sample_info.txt",header=TRUE, row.name=2,sep='\t')
sample_info$title <- factor(sample_info$title,levels=c("WT","KO"))

Enrichment_res <- c("Setoguchi_Tcell_GAGE_KEGG_Biological_function.txt","Setoguchi_Tcell_GAGE_KEGG_Functional_category.txt","Setoguchi_Tcell_GAGE_KEGG_Pathway.txt","Setoguchi_Tcell_GAGE_KEGG_KO.txt")

for (er in Enrichment_res){

	df <- read.table(er,header=TRUE, row.name=1,sep='\t')

	pval <- c()
	score <- c()

	for (i in c(1:nrow(df))){
		data <- na.omit(data.frame(t(df[i,]),sample_info[rownames(t(df[i,])),]))
		colnames(data) <- c("value",colnames(sample_info))
		e <- try(est <- lm( value ~ title,data=data ), silent=FALSE)
		score <- c(score,coef(summary(est))[2,1])
		pval <- c(pval,coef(summary(est))[2,4])
	}

	df_ext <- data.frame(df,score=score,pval=pval,FDR=p.adjust(pval,method="BH"))

	write.table(df_ext,file=sub(".txt","_ex.txt",er),sep="\t")
}