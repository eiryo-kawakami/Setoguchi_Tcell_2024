library(Rsubread)
library(edgeR)

sample.info = read.table('sample_info.txt',header=TRUE,sep='\t')

sample.list = sample.info$sample_id

cwd <- getwd()

files <- list.files()

for (sampleID in sample.list){
	expID = sample.info[sample.info$sample_id==sampleID,]$exp_id
	ispaired = sample.info[sample.info$sample_id==sampleID,]$is_paired
	outputfile = paste(sampleID,"_rc.txt",sep="")
	rpkm = paste(sampleID,"_rpkm.txt",sep="")
	inputfile = paste(sampleID,".Aligned.sortedByCoord.out.bam",sep="")
	if (!(outputfile %in% files)){
		if (ispaired == 'paired'){
			fc <- featureCounts(files=inputfile, annot.ext="/mnt/qnap/ref/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-current/Genes/genes.gtf", isGTFAnnotationFile=TRUE, GTF.attrType="gene_id", isPairedEnd=TRUE, nthreads=24)
		}else{
			fc <- featureCounts(files=inputfile, annot.ext="/mnt/qnap/ref/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-current/Genes/genes.gtf", isGTFAnnotationFile=TRUE, GTF.attrType="gene_name", isPairedEnd=FALSE, nthreads=24)
		}
		x <- DGEList(counts=fc$counts,genes=fc$annotation[,c("GeneID","Length")])
		x_rpkm <- rpkm(x,x$genes$Length)
		write.table(fc$counts, outputfile, quote=F, sep="\t")
		write.table(x_rpkm, rpkm, quote=F, sep="\t")
	}
}
