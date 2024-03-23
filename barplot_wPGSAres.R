library(ggplot2)
library(reshape2)

wPGSAres <- read.table('Setoguchi_Tcell_logFC_TF_wPGSA_t_score.txt',sep="\t",header=T)
top_res <- wPGSAres[c(head(order(wPGSAres$log2FoldChange),5),head(rev(order(wPGSAres$log2FoldChange)),5)),]
top_res <- top_res[order(top_res$log2FoldChange),]
top_res <- top_res[,c(1,5)]
top_res$class <- factor(c(rep("down",5),rep("up",5)),levels=c("up","down"))
top_res$TF <- factor(top_res$TF,levels=top_res$TF)

cols = c("#C7243A","#3261AB")

ggplot(data=top_res)+
geom_bar(aes(x=TF,y=log2FoldChange,fill=class),stat = "identity")+
scale_fill_manual(values=cols)+
theme_bw(base_size = 16) +
theme(
  panel.background = element_rect(size=0.5,colour="black",fill="white")
)+
coord_flip()
ggsave("Setoguchi_Tcell_logFC_TF_wPGSA_t_score_topTF.pdf",useDingbats=F,width=5,height=4)
