# To prepare libraries for this pipeline:
#   1. poretools fasta --type 2D <path/to/base/called/reads/> > <name.fasta>
#   2. bwa mem -x on2d <indexed_reference.fasta> <name.fasta> | samtools view -bS - | samtools sort -o <name.sorted.bam> -
#   3. samtools depth <name.sorted.bam> > <name.coverage>
#   3a.head <name.coverage> # This finds the name of the 'chromosome'; there may be >1.
#   4. awk '$1 == "<chromosomename>" {print $0}' <name.coverage> > chr1.coverage
#   5. Repeat for paired library
#   6. Fill in <name1> and <name2> into pool1 and pool2 below
#
#   Pairs NB01,NB07 Done
#         NB02,NB08 Done
#         NB03,NB09
#         NB04,NB10
#         NB05,NB11
#         NB06,NB12

roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

makeOverlapGraphs <- function(pool1,pool2) {
  
  path1 <- paste('/Users/bpotter/zika-seq-local/',pool1,'/processed/chr1.coverage',sep='')
  path2 <- paste('/Users/bpotter/zika-seq-local/',pool2,'/processed/chr1.coverage',sep='')
  pngName <- paste('/Users/bpotter/zika-seq/depth_coverage/figures/Coverage-Overlap-',pool1,'-',pool2,'.png',sep='')
  
  n1.chr1 <- read.table(path1, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  n2.chr1 <- read.table(path2, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  
  library(reshape)
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
  
  print(round(p20,digits = 4))
  print(round(p40,digits = 4))
  
  png(file=pngName,width=1200,height=600)
  plot(x=n1.chr1$locus, y=n1.chr1$depth, type='l', xlab='locus', ylab='depth', main="Depth of Coverage - Pass reads only", ylim=c(0, graphHeight))
  lines(x=n2.chr1$locus,y=n2.chr1$depth,col="green")
  abline(a=20,b=0,col="orange",lwd=0.5)
  abline(a=20,b=0,col="red",lwd=0.75)
  grid()
  legend(0,max(n1.chr1$depth),c(pool1,pool2),lty=c(1,1),lwd=c(2.5,2.5),col=c('black','green'))
  dev.off()
  rm(list=ls())
}

makeOverlapGraphs('NB01','NB07')
makeOverlapGraphs('NB02','NB08')
makeOverlapGraphs('NB03','NB09')
makeOverlapGraphs('NB04','NB10')
makeOverlapGraphs('NB05','NB11')
makeOverlapGraphs('NB06','NB12')



