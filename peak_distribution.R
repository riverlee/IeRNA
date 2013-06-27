#!/usr/bin/Rscript
##############################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Tue Jun 25 11:33:32 2013
##############################################
# Take a peak file in bed/broadpeak format, and draw width distribution, peaks numbers in each chromosome.
args<-commandArgs(trailingOnly=TRUE)

if(!require("gtools")){
  cat("Package 'gtools' is not found \n")
  q(status=1)
}
#############################################
# Some functions
#############################################
usage<-function(){
    cat("Usage: peak_distribution.R in.bed out.pdf\n");
    q(status=1);
}


if(length(args)!=2){
    usage();
}

in.bed <- args[1]
out.pdf<- args[2]

# Check whether file exists
if(!file.exists(in.bed)){
    cat(paste("File '",in.bed,"' doesn't exists\n",sep=""))
    q(status=1)
}

# read the bed format file
bed.dat<-read.table(in.bed,sep="\t",stringsAsFactors=FALSE,header=FALSE)
cat("Total Number of peaks",nrow(bed.dat),"\n")

# 1) Number of peaks on each chromosome
number.peak.each.chr<-table(bed.dat[,1])
number.peak.each.chr<-number.peak.each.chr[mixedorder(names(number.peak.each.chr))]
cat("Number of peaks on each chromosome\n")
print(number.peak.each.chr)

# 2) Peaks width distibubtion 
bed.dat$width<-abs(bed.dat[,3]-bed.dat[,2])
cat("peaks width summary\n")
print(summary(bed.dat$width))


# 3) plot figures

pdf(file=out.pdf,width=10)
barplot(number.peak.each.chr,space=0,main="Number of peaks on each chromosome",col=rainbow(length(number.peak.each.chr)))

par(mfrow=c(1,2))
boxplot(bed.dat$width,outline=FALSE,main="Width of Peaks (exclude outline)",ylab="Width(bp)")
boxplot(bed.dat$width,outline=TRUE,main="Width of Peaks (include outline)",ylab="Width(bp)")

# each chromosomes
for(chr in names(number.peak.each.chr)){
  par(mfrow=c(1,2))
  boxplot(bed.dat$width[bed.dat$V1==chr],outline=FALSE,main=paste("Width of Peaks on ",chr,"\n (exclude outline,n=",number.peak.each.chr[chr],")",sep=""),ylab="Width(bp)")
  boxplot(bed.dat$width[bed.dat$V1==chr],outline=TRUE,main=paste("Width of Peaks on ",chr,"\n (include outline,n=",number.peak.each.chr[chr],")",sep=""),ylab="Width(bp)")
}

dev.off()
