samplename=$1

echo "${commondir}/${samplename}_H1_L001_R1.fastq.gz"

REF=/mnt/qnap/ref/Mus_musculus/STAR/UCSC/mm10
commondir="/mnt/qnap/Setoguchi_Tcell/sequence/"

CmmnPrms="--runThreadN 24 --outSJfilterReads Unique --outFilterType BySJout --outSAMunmapped Within \
--outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --alignIntronMin 20 \
--outFilterMismatchNoverReadLmax 0.04 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1"
AdtlPrms="--outSAMtype BAM SortedByCoordinate --outBAMcompression 10 --limitBAMsortRAM 32000000000 \
--quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeBAMcompression 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx"

STAR $CmmnPrms $AdtlPrms --genomeDir $REF --outFileNamePrefix ${samplename}. --readFilesIn "${commondir}/${samplename}_H1_L001_R1.fastq.gz" "${commondir}/${samplename}_H1_L001_R2.fastq.gz" --readFilesCommand zcat
