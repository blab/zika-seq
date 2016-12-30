### Overlap Graphs
#### To prepare libraries for this pipeline:
1. poretools fasta --type 2D <path/to/base/called/reads/> > <name.fasta>
2. bwa mem -x on2d <indexed_reference.fasta> <name.fasta> | samtools view -bS - | samtools sort -o <name.sorted.bam> -
3. samtools depth <name.sorted.bam> > <name.coverage>
  - head <name.coverage> # This finds the name of the 'chromosome'; there may be >1.
4. awk '$1 == "<chromosomename>" {print $0}' <name.coverage> > chr1.coverage
5. Repeat for paired library
6. Fill in <name1> and <name2> into pool1 and pool2 below

#### Figures
##### NB01-NB07 Overlap
![](figures/Coverage-Overlap-NB01-NB07.png)

##### NB01-NB07 Overlap
![](figures/Coverage-Overlap-NB02-NB08.png)

##### NB01-NB07 Overlap
![](figures/Coverage-Overlap-NB03-NB09.png)

##### NB01-NB07 Overlap
![](figures/Coverage-Overlap-NB04-NB10.png)

##### NB01-NB07 Overlap
![](figures/Coverage-Overlap-NB05-NB11.png)

##### NB01-NB07 Overlap
![](figures/Coverage-Overlap-NB06-NB12.png)
