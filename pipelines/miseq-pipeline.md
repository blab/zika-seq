# Pipeline for generating consensus sequences from MiSeq run

## Data Sync

Data from a completed MiSeq run will be transferred from the Genomics Core to the lab directory on rhino. Data is then transferred to `/fh/fast/bedford_t/zika-seq/data/` to a library specific directory.

Data should also be copied onto the sequencing data external harddrive `Volumes/Meristem/data/<library>`

## Bioinformatic pipeline

Here we are using a slightly modified version of the [Andersen lab's bioinformatic pipeline](https://github.com/andersen-lab/zika-pipeline).

Install the following (easiest install is with [Homebrew](http://brew.sh/))

    brew install samtools
    brew install trimmomatic
    brew install novoalign

You can have a look at the arguments with `man samtools`, `trimmomatic --help`, and `novoalign --help`.

### Trim off primer sequences

Rather than removing primer sequences based on string matching with the actual primer sequences we trim the first 22 bases from the 5' ends of the reads.

We also trim portions of the leading and trailing ends of the reads if those areas do not align well. Leading and trailing regions are trimmed until the quality score of the read exceeds 20.

Command:

`trimmomatic PE -threads 16  {input.read1} {input.read2} {output.trim1} {output.trim_unpaired1} {output.trim2} {output.trim_unpaired2} ILLUMINACLIP:{input.primer}:2:30:12 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:30 HEADCROP:22`

Example with sample number 569, pool 1:

`trimmomatic PE 569-1_S12_L001_R1_001.fastq.gz 569-1_S12_L001_R2_001.fastq.gz 569-1_R1.trimmed.fastq 569-1_R1.trim_unpaired1.fastq 569-1_R2.trimmed.fastq 569-1_R2.trim_unpaired2.fastq ILLUMINACLIP:amplicons_prefix.fa:2:30:12 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:30 HEADCROP:22`

The vast majority of reads should be in your trimmed fastq files, however because some reads will be unpaired you need to specify an output file for these reads, but they will not be used downstream in the pipeline.

* PE specifies paired-end mode.
* ILLUMINACLIP removes adaptor sequence using palindrome clipping. From the [trimmomatic documentation](http://www.usadellab.org/cms/index.php?page=trimmomatic):

> 'Palindrome' trimming is specifically designed for the case of 'reading through' a short fragment into the adapter sequence on the other end. In this approach, the appropriate adapter sequences are 'in silico ligated' onto the start of the reads, and the combined adapter+read sequences, forward and reverse are aligned. If they align in a manner which indicates 'read-through', the forward read is clipped and the reverse read dropped (since it contains no new data).

* LEADING and TRAILING species trimming the ends of reads until read quality score reaches 20 or higher.
* SLIDINGWINDOW specifies the number of bases to average across, here 4 bases and a required quality score of 25.
* MINLEN will drop any reads that are shorter than the specified integer post crop
* HEADCROP chomps the first 22 bases, this is what removes the primers.
_You can also specify the number of threads to use with_ `-threads`.

**Trimming occurs in the order of the arguments in the command.**

### Map reads to reference

Align reads to reference ZIKV genome Zika (GenBank ID: KU853012) using novoalign. Then pipe the output to samtools to get some basic stats about the reads.

Command:

`novoalign -f {input[0]} {input[1]} -c 16 -r Random -l 40 -g 40 -x 20 -t 502 -d {input[2]} -o SAM | samtools view -F 4 -Sb -o {output}`

Example using sample 569, pool 1, read 1 and read 2:

`novoalign -f 569-1_R1.trimmed.fastq 569-1_R2.trimmed.fastq -r Random -l 40 -g 40 -x 20 -t 502 -d zika_dc_2016.nix -o SAM | samtools view -F 4 -Sb -o 569-1.trimmed.aligned.bam`

_Note: if you're not using samtools you could output to BAM format here instead._

* `-f` Files containing the read sequences to be aligned. By specifying two files the program treats them as paired end reads.
* `-c` specifies the number of threads to use
* `-r` sets the rules for handling of reads with multiple alignment locations. Here `Random` indicates that a single alignment location is randomly chosen from amongst all possible alignment locations.
* `-l` sets the minimum read length.
* `-g` sets gap opening penalty (we are using default of 40)
* `-x` sets gap extension penalty (default is 6, so we are using a higher penalty)
* `-t` sets the maximum alignment score acceptable for the best alignment.
* `-o` specifies output file type
