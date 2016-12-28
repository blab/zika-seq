# Experimental protocols and bioinformatic pipelines for Zika genome sequencing

## Install

Clone the repo:

    git clone https://github.com/blab/zika-seq.git

Install Python dependencies:

    pip install -r requirements.txt

## Data sync

From `zika-seq` run:

    rsync -azP tbedford@rhino.fhcrc.org:/fh/fast/bedford_t/zika-seq/data/ data/

Replacing `tbedford` with your username.

## Bioinformatic pipeline

Data lives in the [`data/`](data/) directory and is not versioned within the repo. Directory structure described in its [README.md](data/).

### Base calling

Convert raw MinION output to FAST5

    metrichor-cli -a <API KEY> -w 1289 -f - -i <directory_with_fast5s> -o downloads

### Convert to FASTA/FASTQ

Extract basecalled information from nanopore FAST5 files

    poretools fasta --type 2D pass/
    poretools fasta --type 2D --highquality fail/
