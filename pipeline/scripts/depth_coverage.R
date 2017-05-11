#!/usr/bin/env Rscript

# To prepare libraries for this pipeline:
#   1. poretools fasta --type 2D <path/to/base/called/reads/> > <name.fasta>
#   2. bwa mem -x on2d <indexed_reference.fasta> <name.fasta> | samtools view -bS - | samtools sort -o <name.sorted.bam> -
#   3. samtools depth <name.sorted.bam> > <name.coverage>
#   3a.head <name.coverage> # This finds the name of the 'chromosome'; there may be >1
#   4. awk '$1 == "<chromosomename>" {print $0}' <name.coverage> > chr1.coverage
#

library(optparse)
library(reshape)

Sys.setenv("DISPLAY"=":0.0")

option_list <- list(
  make_option('--inFile', type='character', help='path to <name.chr1.coverage>'),
  make_option('--outPath', type='character', help='path to the output directory for .png files'),
  make_option('--name', type='character', help='run name')
  )

 parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)

roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

makeOverlapGraphs <- function(infile,odir,runName) {

  path1 <- paste(infile)
  print(infile)
  pngName <- paste(odir,runName,'-Coverage-Overlap.png',sep='')
  logFile <- paste(odir,runName,'-log.txt',sep='')

  n1.chr1 <- read.table(path1, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

  n1.chr1<-rename(n1.chr1,c(V1="Chr", V2="locus", V3="depth"))

  maxDepth <- max(n1.chr1$depth)
  graphHeight <- roundUpNice(maxDepth * 1.1)
  graphHeight

  p20 <- 0
  p40 <- 0
  shorter = max(n1.chr1$locus)
  for (i in 1:shorter) {
    a = n1.chr1$depth[i]
    if (!is.na(a)) {
      if (a >= 40) {
        p20 <- p20 + (1/shorter)
        p40 <- p40 + (1/shorter)
      } else if (a >= 20) {
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
  plot(x=n1.chr1$locus, y=n1.chr1$depth, type='l', xlab='locus', ylab='depth', col='#781C86', main="Depth of Coverage - Pass reads only", ylim=c(0, graphHeight))
  abline(a=40,b=0,col="#D3AE4E",lwd=0.5)
  abline(a=20,b=0,col="#DF4327",lwd=0.75)
  grid()
  legend(0,(.95 * graphHeight),c('Depth','40 read depth','20 read depth'),lty=c(1,1,1),lwd=c(2.5,1,1.5),col=c('#781C86','#D3AE4E','#DF4327'))
  graphics.off()
  rm(list=ls())
}

main <- function () {
  args <- parse_args(parser)

  i <- args$inFile
  print(i)
  print(typeof(i))
  o <- args$outPath
  n <- args$name

  makeOverlapGraphs(i,o,n)
}

main()
