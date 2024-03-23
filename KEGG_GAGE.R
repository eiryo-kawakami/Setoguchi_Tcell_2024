library(ggplot2)
library(gage)



logfc.file = "Setoguchi_Tcell_logFC_gID.txt"

KEGG_category = c("KO","Pathway","Biological_function","Functional_category")

for (cat in KEGG_category){

	exp.fc <- read.table(logfc.file,sep='\t',header=T,row.names=1)
	exp.fc <- exp.fc[,c(2),drop=F]

	gset_file <- readList(paste("KEGG_",cat,"_gset.gmt",sep=""))

	kegg.p <- gage(exp.fc,gsets=gset_file, set.size = c(10, 1000))
	kegg.sig <- sigGeneSet(kegg.p,outname=paste("Setoguchi_Tcell_KEGG_",cat,sep=""),cutoff=1)
	write.table(rbind(kegg.sig$greater,kegg.sig$less),file=paste("Setoguchi_Tcell_GAGE_KEGG_",cat,".txt",sep=""),sep='\t',quote=T)

}
