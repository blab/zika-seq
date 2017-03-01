# To prepare libraries for this pipeline:
#   1. poretools fasta --type 2D <path/to/base/called/reads/> > <name.fasta>
#   2. bwa mem -x on2d <indexed_reference.fasta> <name.fasta> | samtools view -bS - | samtools sort -o <name.sorted.bam> -
#   3. samtools depth <name.sorted.bam> > <name.coverage>
#   3a.head <name.coverage> # This finds the name of the 'chromosome'; there may be >1
#   4. awk '$1 == "<chromosomename>" {print $0}' <name.coverage> > chr1.coverage
#   5. Repeat for paired library
#   6. Fill in <name1> and <name2> into pool1 and pool2 below
#
#   Pairs NB01,NB07
#         NB02,NB08
#         NB03,NB09
#         NB04,NB10
#         NB05,NB11
#         NB06,NB12

library(argparse)
library(reshape)

parser <- ArgumentParser(description='Make coverage graphs from multiplexed minION reads.')
parser$add_argument('-i','--inputPath', type='character', help='path to directory contining multiplexed reads')
parser$add_argument('-o','--outputPath', type='character', help='path to the output directory for .png files')

roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

makeOverlapGraphs <- function(pool1,pool2,idir,odir) {
  
  path1 <- paste(idir,pool1,'/processed/chr1.coverage',sep='')
  path2 <- paste(idir,pool2,'/processed/chr1.coverage',sep='')
  pngName <- paste(odir,'Coverage-Overlap-',pool1,'-',pool2,'.png',sep='')
  logFile <- paste(odir,'log',pool1,'-',pool2,'.txt',sep='')
  
  n1.chr1 <- read.table(path1, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  n2.chr1 <- read.table(path2, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  
  n1.chr1<-rename(n1.chr1,c(V1="Chr", V2="locus", V3="depth"))
  n2.chr1<-rename(n2.chr1,c(V1="Chr", V2="locus", V3="depth"))
  
  maxDepth <- max(max(n1.chr1$depth),max(n2.chr1$depth))
  graphHeight <- roundUpNice(maxDepth * 1.1)
  graphHeight
  
  p20 <- 0
  p40 <- 0
  shorter = min(max(n1.chr1$locus),max(n2.chr1$locus))
  for (i in 1:shorter) {
    a = n1.chr1$depth[i]
    b = n2.chr1$depth[i]
    if (!is.na(a) & !is.na(b)) {
      if ((a + b) >= 40) {
        p20 <- p20 + (1/shorter)
        p40 <- p40 + (1/shorter)
      } else if ((a + b) >= 20) {
        p20 <- p20 + (1/shorter)
      }
    }
  }
  
  p20 <- round(p20,digits = 4)
  p40 <- round(p40,digits = 4)
  write('Percentage over 20x depth:',logFile,append=TRUE)
  write(p20,logFile,append=TRUE)
  write('Percentage over 40x depth:',logFile,append=TRUE)
  write(p40,logFile,append=TRUE)
  
  png(file=pngName,width=1200,height=600)
  plot(x=n1.chr1$locus, y=n1.chr1$depth, type='l', xlab='locus', ylab='depth', main="Depth of Coverage - Pass reads only", ylim=c(0, graphHeight))
  lines(x=n2.chr1$locus,y=n2.chr1$depth,col="yellowgreen")
  abline(a=40,b=0,col="goldenrod2",lwd=0.5)
  abline(a=20,b=0,col="tomato2",lwd=0.75)
  grid()
  legend(0,(.95 * graphHeight),c(pool1,pool2,'40 read depth','20 read depth'),lty=c(1,1,1,1),lwd=c(2.5,2.5,1,1.5),col=c('black','yellowgreen','goldenrod2','tomato2'))
  # dev.off()
  graphics.off()
  rm(list=ls())
}

main <- function () {
  args <- parser$parse_args()
  
  makeOverlapGraphs('NB01','NB07',args$inputPath,args$outputPath)
  makeOverlapGraphs('NB02','NB08',args$inputPath,args$outputPath)
  makeOverlapGraphs('NB03','NB09',args$inputPath,args$outputPath)
  makeOverlapGraphs('NB04','NB10',args$inputPath,args$outputPath)
  makeOverlapGraphs('NB05','NB11',args$inputPath,args$outputPath)
  makeOverlapGraphs('NB06','NB12',args$inputPath,args$outputPath)
}

main()