# Experimental protocols and bioinformatic pipelines for Zika genome sequencing

## Install

Clone the repo:

    git clone https://github.com/blab/zika-seq.git

Install Python dependencies:

    pip install -r requirements.txt

Install [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library):

    git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.git
    cd Complete-Striped-Smith-Waterman-Library/src/
    cp libssw.so <path to zika-seq>/zika-seq/scripts/
    cp ssw_lib.py <path to zika-seq>/zika-seq/scripts/

Install [marginAlign](https://github.com/benedictpaten/marginAlign):

    git clone https://github.com/benedictpaten/marginAlign.git
    cd marginAlign
    git submodule update --init --recursive
    make
    export PATH=<path to marginAlign>/marginAlign/:$PATH

Install [samtools](https://github.com/samtools/samtools):

    brew tap homebrew/science
    brew install samtools

## Data sync

From `zika-seq` run:

    rsync -azP tbedford@rhino.fhcrc.org:/fh/fast/bedford_t/zika-seq/data/ data/

Replacing `tbedford` with your username.

## Bioinformatic pipeline

Data lives in the [`data/`](data/) directory and is not versioned within the repo. Directory structure described in its [README.md](data/).

### Base calling

Convert raw MinION output to FAST5

    metrichor-cli -a <API KEY> -w 1289 -f - -i <directory_with_fast5s> -o downloads

### Run pipeline

Run poretools, marginAlign, samtools:

    python run.py
