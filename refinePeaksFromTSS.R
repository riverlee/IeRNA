#!/usr/bin/Rscript
##############################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Tue Jun 25 13:37:40 2013
##############################################
args<-commandArgs(trailingOnly=TRUE)

#############################################
# Some functions
#############################################
usage<-function(){
  cat("Usage: refinePeaksFromTSS.R peak.bed refinedTSS.bed 1000 refinedPeaks.bed\n");
  q(status=1);
}

if(length(args)!=4){
  usage();
}

in.bed <- args[1]
tss.bed <- args[2]
dist<- as.numeric(args[3])
out.bed<- args[4]

# Check whether file exists
if(!file.exists(in.bed)){
  cat(paste("File '",in.bed,"' doesn't exists\n",sep=""))
  q(status=1)
}

# Check whether file exists
if(!file.exists(tss.bed)){
  cat(paste("File '",tss.bed,"' doesn't exists\n",sep=""))
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

# 1) Read the peak bed file
peak.dat<-read.table(in.bed,sep="\t",stringsAsFactors=FALSE,header=FALSE)
cat("Total Number of peaks",nrow(peak.dat),"\n")
peak.range<-GRanges(seqnames=peak.dat[,1],ranges=IRanges(start=peak.dat[,2],end=peak.dat[,3]))

# 2) Read TSS 
tss.dat<-read.table(tss.bed,sep="\t",stringsAsFactors=FALSE,header=FALSE)
cat("Total Number of TSSs",nrow(tss.dat),"\n")
tss.range<-GRanges(seqnames=tss.dat[,1],ranges=IRanges(start=tss.dat[,2],end=tss.dat[,3]))
# Base on strand, get the 5' position
tss.start<-apply(X=tss.dat,MARGIN=1,function(x){ifelse(x[6]=="+",x[2],x[3])})
tss.start<-as.numeric(tss.start)
tss.range2<-GRanges(seqnames=tss.dat[,1],ranges=IRanges(start=tss.start,end=tss.start+1))

# 3) Exclude those peaks overlapped with gene body
cat("Exclude those peaks overlapped with gene body ....\n")
ol<-findOverlaps(query=peak.range,subject=tss.range)
ol<-as.matrix(ol)
cat("Exclude ",length(unique(ol[,1])),"Peaks\n")

cat("Exclude those peaks within ",dist,"bp of gene's TSS ....\n",sep="")
ol2<-findOverlaps(query=peak.range,subject=tss.range2,maxgap=dist)
ol2<-as.matrix(ol2)
ol2<-setdiff(unique(ol2[,1]),unique(ol[,1]))
cat("Exclude another",length(unique(ol2)),"Peaks\n")

# Get the index of query
peak.dat.remained<-peak.dat[-c(unique(ol[,1]),ol2),]
cat("Remaining ",nrow(peak.dat.remained),"Peaks\n")

# Write out
write.table(x=peak.dat.remained,file=out.bed,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

