#Rscript pathview_new.R Stat3KO 10 5 contrast.txt など

library(pathview)
library(dplyr)

logfc.file = "Setoguchi_Tcell_logFC_gID.txt"
exp.fc <- read.table(logfc.file,sep='\t',header=T,row.names=1)

tbl <- read.table("Setoguchi_Tcell_GAGE_KEGG_Pathway_ex.txt", header=T, row.names=1)
dec <- tbl[tbl$pval<0.01,]
path.ids <- c(rownames(dec))
path.ids <- apply(path.ids,c(1,2),sub(pattern="map",replacement="mmu"))
path.ids <- path.ids %>% sub(pattern="map",replacement="mmu")

for (pid in path.ids){
	pv.out.list <- pathview(gene.data = exp.fc[,2,drop=F],pathway.id = substr(pid,1,8), species="mmu",out.suffix=substr(gsub(" ","_",pid),10,100))
	file.remove(paste(substr(pid,1,8),".png",sep=""))
	file.remove(paste(substr(pid,1,8),".xml",sep=""))
}
