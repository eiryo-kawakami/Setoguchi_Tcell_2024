library(annotate)
library(org.Mm.eg.db)

genesym_files <- c("Setoguchi_Tcell_logFC.txt","Setoguchi_Tcell_rlog_count_center_cutoff.txt")

for (f in genesym_files){

	geneID_file <- sub(".txt$","_gID.txt",f)
	x <- read.table(f,header=T,row.names=1)
	Gene_Sym <- rownames(x)
	fc <- x[,1:ncol(x)]
	Gene_Sym <- as.character(Gene_Sym)
	GeneID <- select(org.Mm.eg.db,Gene_Sym,"ENTREZID","SYMBOL")
	GeneID <- GeneID[!duplicated(GeneID$SYMBOL),]$ENTREZID
	y <- data.frame(Gene_Sym=Gene_Sym,GeneID=GeneID,fc)
	y <- y[!is.na(y$GeneID),]
	write.table(y,geneID_file,quote=F,col.names=T,row.names=F,sep="\t")

}