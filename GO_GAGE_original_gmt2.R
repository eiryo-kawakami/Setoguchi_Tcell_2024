# Rscript KEGG_GO_GAGE_new.R contrast.txt Stat3KO

library(ggplot2)
library(gage)
#data(kegg.sets.mm)

logfc.file = "Setoguchi_Tcell_logFC_gID.txt"

GO_category = c("cellular_component","biological_process","molecular_function")

for (cat in GO_category){

	exp.fc <- read.table(logfc.file,sep='\t',header=T,row.names=1)
	exp.fc <- exp.fc[,c(2),drop=F]

	gset_file <- readList(paste("GO_gset_",cat,".gmt",sep=""))

	go.p <- gage(exp.fc,gsets=gset_file, set.size = c(10, 1000))
	go.sig <- sigGeneSet(go.p,outname=paste("Setoguchi_Tcell_GO_",cat,sep=""),cutoff=1)
	write.table(rbind(go.sig$greater,go.sig$less),file=paste("Setoguchi_Tcell_GAGE_GO_",cat,".txt",sep=""),sep='\t',quote=T)

}
