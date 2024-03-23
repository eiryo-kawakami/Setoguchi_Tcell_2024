library(DESeq2)
library(ggplot2)
library(genefilter)
library(Rtsne)
library(data.table)
library(RColorBrewer)
library(cluster)
library(ggrepel)
library(matrixStats)

limit <- 10

cutoff <- function(num,limit){
	if(num > limit){
		return(limit)
	}else if(num < -limit){
		return(-limit)
	}else return(num)
}

if (exists("pamNew") ) rm(pamNew)
pamNew <- function (x, k, diss1 = inherits(x, "dist"), metric1 = "euclidean")
{
  
  #############################################################################################################
  # A new pam clustering function which corrects the clustering membership based on the sillhouette strength. #
  # The clustering membership of an observation with a negative sillhouette strength is reassigned to its     #
  # neighboring cluster.                                                                                      #
  # The inputs of the function are similar to the original 'pam' function.                                    #
  # The function returns a vector of clustering labels.                                                       #
  # Copyright 2003 Tao Shi and Steve Horvath (last modified 10/31/03)                                         #
  #############################################################################################################
  
  if (diss1)
  {
    if (!is.null(attr(x, "Labels"))) { original.row.names <- attr(x, "Labels")}
    names(x) <- as.character(c(1:attr(x, "Size")))
  } 
  else
  {
    if(!is.null(dimnames(x)[[1]])) { original.row.names <- dimnames(x)[[1]]}
    row.names(x) <- as.character(c(1:dim(x)[[1]]))
  }
  pam1 <- pam(x,k,diss=diss1, metric=metric1)
  label2 <- pam1$clustering
  silinfo1 <- pam1$silinfo$widths
  index1 <- as.numeric(as.character(row.names(silinfo1)))
  silinfo2 <- silinfo1[order(index1),]
  labelnew <- ifelse(silinfo2[,3]<0, silinfo2[,2], silinfo2[,1])
  names(labelnew) <- original.row.names
  labelnew    
}


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
# mat <- mat - rowMedians(as.matrix(mat))

# mat <- apply(mat,c(1,2),cutoff,limit)

#mat <- mat[,as.vector(rownames(sample_info[sample_info$disease!="healthy",]))]

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
#scale_shape_manual(values = shapes) +
#stat_ellipse(aes(x=tSNE1,y=tSNE2,fill=cluster),geom="polygon", level=0.95, alpha=0.2) +
#geom_path(aes(linetype=pt_ID,color=pt_ID),arrow = arrow(angle = 15,length=unit(0.30,"cm"),type = "closed")) +
#scale_linetype_manual(values=linetypes)+
# scale_x_continuous(limits = c(-200,300)) +
# scale_y_continuous(limits = c(-100,200)) +
theme_bw(base_size = 16) +
theme(
  panel.background = element_rect(size=0.5,colour="black",fill="white"),
)
ggsave("Setoguchi_Tcell_PCA.pdf",useDingbats=FALSE,width=5.7,height=4)


# rlog_mat_scale <- t(mat)

# logFC <- read.table('Ishikawa_AML_disease_AML_vs_healthy_logFC.txt',sep="\t")
# gene_list <- rownames(logFC[head(order(logFC$pvalue),100),])

# for (gene in gene_list){

# 	pca.data.mat <- data.frame(pca.data,gene_exp=rlog_mat_scale[rownames(pca.data),gene])
# 	p <- ggplot(aes(x=PC1,y=PC2),data=pca.data.mat) +
# 	geom_point(aes(colour=gene_exp),shape=16,alpha=0.8,size=3) +
# 	#scale_colour_gradientn(colours=rainbow(4)) + 
# 	scale_colour_distiller(palette = "Spectral") +
# 	#scale_shape_manual(values=c(1,17))+
# 	#stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
# 	#              geom="polygon", level=0.95, alpha=0.2) +
# 	#scale_fill_manual(values=rev(cols[c(1,2)]))+
# 	#scale_colour_manual(values=rev(cols[c(1,2)])) +
# 	scale_x_continuous(limits = c(-200,300)) +
# 	scale_y_continuous(limits = c(-100,200)) +
# 	theme_bw(base_size = 16) +
# 	theme(
# 		panel.background = element_rect(size=0.5,colour="black",fill="white")
# 	)
# 	ggsave(file = paste("Ishikawa_AML_PCA_",gene,".pdf",sep=""), plot = p,useDingbats=FALSE,width=5,height=4)
# }


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

# p <- ggplot(aes(x=tSNE1,y=tSNE2),data=tsne.data) +
# geom_point(aes(color=as.factor(cluster)),shape=16,alpha=0.8,size=3) +
# stat_ellipse(aes(x=tSNE1,y=tSNE2,color=as.factor(cluster)), level=0.95, linetype=2) +
# #scale_colour_brewer(palette = "Set1") + 
# #scale_shape_manual(values=c(1,17))+
# #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
# #              geom="polygon", level=0.95, alpha=0.2) +
# #scale_fill_manual(values=rev(cols[c(1,2)]))+
# scale_fill_manual(values=cols) +
# scale_colour_manual(values=cols) +
# # scale_x_continuous(limits = c(-60,45)) +
# # scale_y_continuous(limits = c(-60,60)) +
# theme_bw(base_size = 16) +
# theme(
#   panel.background = element_rect(size=0.5,colour="black",fill="white"),
# )
# ggsave(file = "NIID_Kato_Virus_tSNE_cluster.pdf", plot = p,useDingbats=FALSE)

ggplot(aes(x=tSNE1,y=tSNE2),data=tsne.data) +
geom_point(aes(color=title),size=3,alpha=0.8,shape=16) +
scale_color_manual(values=cols) +
#scale_shape_manual(values = shapes) +
#stat_ellipse(aes(x=tSNE1,y=tSNE2,fill=cluster),geom="polygon", level=0.95, alpha=0.2) +
#geom_path(aes(linetype=pt_ID,color=pt_ID),arrow = arrow(angle = 15,length=unit(0.30,"cm"),type = "closed")) +
#scale_linetype_manual(values=linetypes)+
# scale_x_continuous(limits = c(-60,45)) +
# scale_y_continuous(limits = c(-60,60)) +
theme_bw(base_size = 16) +
theme(
  panel.background = element_rect(size=0.5,colour="black",fill="white"),
)
ggsave("Setoguchi_Tcell_tSNE.pdf",useDingbats=FALSE,width=7,height=4)

# rlog_mat_scale <- mat

# logFC <- read.table('Ishikawa_AML_disease_AML_vs_healthy_logFC.txt',sep="\t")
# gene_list <- rownames(logFC[head(order(logFC$pvalue),100),])

# for (gene in gene_list){

# 	tsne.data.mat <- data.frame(tsne.data,gene_exp=rlog_mat_scale[rownames(tsne.data),gene])
# 	p <- ggplot(aes(x=tSNE1,y=tSNE2),data=tsne.data.mat) +
# 	geom_point(aes(colour=gene_exp),shape=16,alpha=0.8,size=3) +
# 	#scale_colour_gradientn(colours=rainbow(4)) + 
# 	scale_colour_distiller(palette = "Spectral") +
# 	#scale_shape_manual(values=c(1,17))+
# 	#stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
# 	#              geom="polygon", level=0.95, alpha=0.2) +
# 	#scale_fill_manual(values=rev(cols[c(1,2)]))+
# 	#scale_colour_manual(values=rev(cols[c(1,2)])) +
# 	# scale_x_continuous(limits = c(-60,45)) +
# 	# scale_y_continuous(limits = c(-60,60)) +
# 	theme_bw(base_size = 16) +
# 	theme(
# 		panel.background = element_rect(size=0.5,colour="black",fill="white")
# 	)
# 	ggsave(file = paste("Ishikawa_AML_tSNE_",gene,".pdf",sep=""), plot = p,useDingbats=FALSE,width=7,height=6)
# }
