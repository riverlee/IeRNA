#!/usr/bin/Rscript
##############################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Wed Jun 26 13:52:49 2013
##############################################
args<-commandArgs(trailingOnly=TRUE)

# Take two peaks, keep those peaks from input1 which has peaks from antoher input with 2kb distance

#############################################
# Some functions
#############################################
usage<-function(){
  cat("Usage: overlapOtherPeaks.R peak.bed anotherPeak.bed 2000 output.peak\n");
  q(status=1);
}

if(length(args)!=4){
  usage();
}

in.bed <- args[1]
another.bed <-args[2]
dist<- as.numeric(args[3])
out.bed<- args[4]

# Check whether file exists
if(!file.exists(in.bed)){
  cat(paste("File '",in.bed,"' doesn't exists\n",sep=""))
  q(status=1)
}

if(!file.exists(another.bed)){
  cat(paste("File '",another.bed,"' doesn't exists\n",sep=""))
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

# 2) Read anohter peak bed file
another.dat<-read.table(another.bed,sep="\t",stringsAsFactors=FALSE,header=FALSE)
cat("Total Number of peaks from another dataset",nrow(another.dat),"\n")
another.range<-GRanges(seqnames=another.dat[,1],ranges=IRanges(start=another.dat[,2],end=another.dat[,3]))

# 3) Get those overlaps
cat("Get peaks overlap with another peak dataset....\n")
ol<-findOverlaps(query=peak.range,subject=another.range,maxgap=dist)
ol<-as.matrix(ol)
cat("There are ",length(unique(ol[,1])),"peaks overlap with another peaks dataset\n")

# Get the index of query
peak.dat.remained<-peak.dat[c(unique(ol[,1])),]
cat("Remaining ",nrow(peak.dat.remained),"Peaks\n")

# Write out
write.table(x=peak.dat.remained,file=out.bed,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


