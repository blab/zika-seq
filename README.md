# Experimental protocols and bioinformatic pipelines for Zika genome sequencing

## Install

Clone the repo:

    git clone https://github.com/blab/zika-seq.git

Install Python dependencies:

    pip install -r requirements.txt

## Bioinformatic pipeline

Data lives in the [`data/`](data/) directory and is not versioned within the repo. Directory structure described in its [README.md](data/).

### Base calling

Convert raw MinION output to FAST5

    metrichor-cli

### Convert to FASTA/FASTQ

Extract basecalled information from nanopore FAST5 files

    poretools fasta --type 2D pass/
    poretools fasta --type 2D --highquality fail/
