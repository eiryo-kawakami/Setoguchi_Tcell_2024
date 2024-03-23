library(annotate)
library(org.Mm.eg.db)

rcID_file = "Setoguchi_Tcell_logFC.txt"
rc_geneID <- sub(".txt$","_gID.txt",rcID_file)
x <- read.table(rcID_file,header=T,row.names=1)
Gene_Sym <- rownames(x)
fc <- x[,1:ncol(x)]
Gene_Sym <- as.character(Gene_Sym)
GeneID <- select(org.Mm.eg.db,Gene_Sym,"ENTREZID","SYMBOL")
GeneID <- GeneID[!duplicated(GeneID$SYMBOL),]$ENTREZID
y <- data.frame(Gene_Sym=Gene_Sym,GeneID=GeneID,fc)
y <- y[!is.na(y$GeneID),]
write.table(y,rc_geneID,quote=F,col.names=T,row.names=F,sep="\t")
