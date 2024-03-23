library(DESeq2)

sample_db <- read.table("sample_info.txt",sep="\t",header=T,row.names=2)
sample_db$title <- factor(sample_db$title,levels=c("WT","KO"))

count_summary <- read.table('Setoguchi_Tcell_readcount_summary.txt',sep="\t",header=T,row.names=1)

count_summary <- count_summary[,as.vector(rownames(sample_db))]
count_summary <- count_summary[rowSums(count_summary)!=0,]
count_summary1 <- sapply(count_summary,as.integer)
rownames(count_summary1) <- rownames(count_summary)

dds <- DESeqDataSetFromMatrix(countData=count_summary1, colData=sample_db, design=~title+source)
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# dds <- nbinomLRT(dds, full = ~ cluster, reduced = ~ 1)
dds <- DESeq(dds, betaPrior=FALSE)

res <- lfcShrink(dds, coef="title_KO_vs_WT", type="apeglm")
#res <- results(dds, name="title_KO_vs_WT",independentFiltering=FALSE)

write.table(res,file="Setoguchi_Tcell_logFC.txt",sep="\t",quote=F)
