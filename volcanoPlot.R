library(ggplot2)
library(ggrepel)

logFC_threshold1 = 0.25
logFC_threshold2 = 0.5
FDR_threshold = 0.05

cols = c("#C7243A","#555555","#3261AB")
#cols = c("#C7243A","#DA6272","#555555","#6A8CC7","#3261AB")

logFCdata <- read.table('Setoguchi_Tcell_logFC_gID.txt',header=T)
logFCdata$padj <- ifelse(is.na(logFCdata$padj),1,logFCdata$padj)
logFCdata$gene_cat <- rep("N",nrow(logFCdata))
#logFCdata$gene_cat <- ifelse(logFCdata$padj<FDR_threshold & logFCdata$log2FoldChange>logFC_threshold1,"H",logFCdata$gene_cat)
logFCdata$gene_cat <- ifelse(logFCdata$padj<FDR_threshold & logFCdata$log2FoldChange>logFC_threshold2,"HH",logFCdata$gene_cat)
#logFCdata$gene_cat <- ifelse(logFCdata$padj<FDR_threshold & logFCdata$log2FoldChange< -logFC_threshold1,"L",logFCdata$gene_cat)
logFCdata$gene_cat <- ifelse(logFCdata$padj<FDR_threshold & logFCdata$log2FoldChange< -logFC_threshold2,"LL",logFCdata$gene_cat)

logFCdata$gene_cat <- factor(logFCdata$gene_cat,levels=c("HH","H","N","L","LL"))
logFCdata$padj <- as.numeric(-log10(logFCdata$padj))
logFCdata$log2FoldChange <- as.numeric(logFCdata$log2FoldChange)

logFCdata_selected_H <- logFCdata[logFCdata$gene_cat=="HH",]
write.table(logFCdata_selected_H,file="Setoguchi_Tcell_logFC_gID_upDEG.txt",sep="\t",quote=F,row.names=F)
logFCdata_selected_L <- logFCdata[logFCdata$gene_cat=="LL",]
write.table(logFCdata_selected_L,file="Setoguchi_Tcell_logFC_gID_downDEG.txt",sep="\t",quote=F,row.names=F)

ggplot()+
geom_point(data=logFCdata,aes(x=log2FoldChange,y=padj,color=gene_cat,size=gene_cat),alpha=0.8,shape=16) +
scale_color_manual(values=cols) +
scale_size_manual(values=c(1.5,1,1.5))+
geom_text_repel(data=logFCdata_selected_H,aes(x=log2FoldChange,y=padj,label=GeneName),color="#C7243A",size=3)+
geom_text_repel(data=logFCdata_selected_L,aes(x=log2FoldChange,y=padj,label=GeneName),color="#3261AB",size=3)+
geom_vline(xintercept =c(-logFC_threshold2,logFC_threshold2),linetype="dashed")+
geom_hline(yintercept =-log10(FDR_threshold),linetype="dashed")+
theme_bw(base_size = 16) +
theme(
  panel.background = element_rect(size=0.5,colour="black",fill="white")
)
ggsave("Setoguchi_Tcell_gID_volcanoPlot.pdf",useDingbats=FALSE,width=7,height=4)

