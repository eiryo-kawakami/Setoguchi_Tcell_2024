library(DESeq2)
library(ggplot2)
library(genefilter)
library(Rtsne)
library(data.table)
library(RColorBrewer)
library(cluster)
library(ggrepel)
library(matrixStats)


linetypes <- rep(c(1:12),10)
cols <-c(
  brewer.pal(8, "Set1")[1],
  "#333333",
  brewer.pal(8, "Set1")[2],
  "#555555",
  brewer.pal(9, "Greens")[9],
  brewer.pal(9, "Greens")[6],
  brewer.pal(9, "Greens")[3],
  brewer.pal(8, "Set1")[4],
  brewer.pal(8, "Set1")[5],
  "#777777",
  brewer.pal(8, "Set1")[7],
  "#999999",
  brewer.pal(8, "Set1")[8],
  "#BBBBBB")
cols2 <- rep(c(brewer.pal(12, "Paired"),"#333333"),20)
cols3 <- c(brewer.pal(7, "Set1")[c(1,3,4,5)],"#777777")
shapes <- rep(c(0:10,12:13,15:18),8)

sample_info <- read.table('sample_info.txt',sep="\t",header=TRUE,row.names=2)

mat <- read.table("Setoguchi_Tcell_rlog_count_center_cutoff.txt")

pca.res <- prcomp(t(mat))
pca.data <- data.frame(PC1=pca.res$x[,1],PC2=pca.res$x[,2])
rownames(pca.data) <- rownames(t(mat))
pca.data <- data.frame(pca.data,sample_info[rownames(pca.data),])

write.table(pca.data,file="Setoguchi_Tcell_PCA_data.txt",sep='\t',quote=F)

pca.data$title <- factor(pca.data$title,levels=c("WT","KO"))


ggplot(aes(x=PC1,y=PC2),data=pca.data) +
geom_point(aes(color=title),size=3,alpha=0.8,shape=16) +
geom_text_repel(aes(label=sample_name)) +
scale_color_manual(values=cols) +
theme_bw(base_size = 16) +
theme(
  panel.background = element_rect(size=0.5,colour="black",fill="white"),
)
ggsave("Setoguchi_Tcell_PCA.pdf",useDingbats=FALSE,width=5.7,height=4)

set.seed(1)
tsne.res <- Rtsne(t(mat),theta=0.0,initial_dims=50,perplexity=2)
tsne.data <- data.frame(tSNE1=tsne.res$Y[,1],tSNE2=tsne.res$Y[,2])
rownames(tsne.data) <- rownames(t(mat))
tsne.data <- data.frame(tsne.data,sample_info[rownames(tsne.data),])

tSNEclusterLabel = pamNew(tsne.data[,c(1,2)], 3)

tsne.data <- data.frame(tsne.data,cluster=tSNEclusterLabel)

write.table(tsne.data,file="Setoguchi_Tcell_tSNE_data.txt",sep='\t',quote=F)

tsne.data <- na.omit(read.table('Setoguchi_Tcell_tSNE_data.txt',sep="\t"))

tsne.data$title <- factor(tsne.data$title,levels=c("WT","KO"))

ggplot(aes(x=tSNE1,y=tSNE2),data=tsne.data) +
geom_point(aes(color=title),size=3,alpha=0.8,shape=16) +
scale_color_manual(values=cols) +
theme_bw(base_size = 16) +
theme(
  panel.background = element_rect(size=0.5,colour="black",fill="white"),
)
ggsave("Setoguchi_Tcell_tSNE.pdf",useDingbats=FALSE,width=7,height=4)
