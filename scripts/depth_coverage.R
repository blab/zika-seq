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

pool1 <- 'NB06'
pool2 <- 'NB12'

path1 <- paste('/Users/bpotter/zika-seq/',pool1,'/processed/chr1.coverage',sep='')
path1
path2 <- paste('/Users/bpotter/zika-seq/',pool2,'/processed/chr1.coverage',sep='')
path2
pngName <- paste('/Users/bpotter/zika-seq/Coverage-Overlap-',pool1,'-',pool2,'.png',sep='')
pngName

n1.chr1 <- read.table(path1, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
n1.chr1
n2.chr1 <- read.table(path2, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
n2.chr1

library(reshape)
n1.chr1<-rename(n1.chr1,c(V1="Chr", V2="locus", V3="depth"))
n2.chr1<-rename(n2.chr1,c(V1="Chr", V2="locus", V3="depth"))

png(file=pngName,width=1200,height=600)
plot(x=n1.chr1$locus, y=n1.chr1$depth, type='l', xlab='locus', ylab='depth', main="Depth of Coverage - Pass reads only")
lines(x=n2.chr1$locus,y=n2.chr1$depth,col="green")
legend(0,max(n1.chr1$depth),c(pool1,pool2),lty=c(1,1),lwd=c(2.5,2.5),col=c('black','green'))
# dev.copy(png,pngName)
dev.off()
rm(list=ls())
