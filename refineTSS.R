#!/usr/bin/Rscript
##############################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Tue Jun 25 12:41:34 2013
##############################################
# Refine tss which is download from UCSC table (like refGene)
# For Genes has multiple annotated TSS (transcripts/isoform),only use 5'-most TSSs from these genes
args<-commandArgs(trailingOnly=TRUE)
#############################################
# Some functions
#############################################
usage<-function(){
  cat("Usage: refineTSS.R refGene.txt refinedTSS.bed\n");
  q(status=1);
}

if(length(args)!=2){
  usage();
}

in.bed <- args[1]
out.bed<- args[2]

# Check whether file exists
if(!file.exists(in.bed)){
  cat(paste("File '",in.bed,"' doesn't exists\n",sep=""))
  q(status=1)
}

if(!require(IRanges)){
  cat("Package 'IRanges' is not found \n")
  q(status=1)
}
if(!require(GenomicRanges)){
  cat("Package 'GenomicRanges' is not found \n")
  q(status=1)
}

# read the bed format file
bed.dat<-read.table(in.bed,sep="\t",stringsAsFactors=FALSE,header=TRUE,comment.char="")

# remove those on gl_random
bed.dat<-bed.dat[!grepl(pattern="_",x=bed.dat$chrom),]

# Make GRanges
gene.index<-bed.dat$name2!=""
gene.gr<-GRanges(seqnames=bed.dat$chrom[gene.index],
                 ranges=IRanges(bed.dat$txStart[gene.index],bed.dat$txEnd[gene.index]),
                 strand=bed.dat$strand[gene.index],
                 tx_id=bed.dat$name[gene.index],
                 gene_id=bed.dat$name2[gene.index])
gene.gr.list<-split(gene.gr,bed.dat$name2[gene.index])

# reduce, similar to collapse
gene.gr.list<-reduce(gene.gr.list)

# For each gene, get the tss start, tss end and strand in bed format
gene.gr.list.start.end<-range(gene.gr.list)

# conver to data frame format
gene.tss.start.end<-as.data.frame(gene.gr.list.start.end)
colnames(gene.tss.start.end)<-c("gene","chr","start","end","width","strand")

gene.tss.start.end<-gene.tss.start.end[,c("chr","start","end","gene","width","strand")]

# Write out the bed file
write.table(x=gene.tss.start.end,file=out.bed,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
