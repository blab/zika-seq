# Experimental protocols and bioinformatic pipelines for Zika genome sequencing

#### Allison Black<sup>1,2</sup>, Barney Potter<sup>2</sup>, Nicholas J. Loman<sup>3</sup>, Trevor Bedford<sup>2</sup>

<sup>1</sup>Department of Epidemiology, University of Washington, Seattle, WA, USA, <sup>2</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA, <sup>3</sup>Institute of Microbiology and Infection, University of Birmingham, Birmingham, UK

## Install

Clone the repo and load submodules:

    git clone https://github.com/blab/zika-seq.git
    git submodule update --init --recursive

## Data sync

Primary sequencing data lives on the Rhino FHCRC cluster at:

    /fh/fast/bedford_t/data/

And locally on Meristem drive at:

    /Volumes/Meristem/data/

To sync Meristem to Rhino, run:

    rsync -azP tbedford@rhino.fhcrc.org:/fh/fast/bedford_t/data/ /Volumes/Meristem/data/

Replacing `tbedford` with your username.

This `data/` directory is assumed to follow [a particular schema](https://github.com/blab/zika-seq/blob/master/data-schema.md).

## Bioinformatic pipeline

Here, we use the ZiBRA project bioinformatic pipeline at [zibraproject/zika-pipeline](https://github.com/zibraproject/zika-pipeline/). This pipeline is instantiated in the Docker image [zibra/zibra](https://hub.docker.com/r/zibra/zibra/). Data processing is done using Docker.

### Data volume

Create a named data volume that mirrors local `data/` to `data/` within container:

    docker create --name zibra-data -v /Volumes/Meristem/data:/data zibra/zibra    

This is to get data into the Docker container. Note that the path to local directory has to be an absolute path.

Create a named data volume for a single sample:

    docker create --name zibra-data-lb01-nb01 -v /Volumes/Meristem/data/usvi-library1-2016-12-10/basecalled_reads/pass_demultiplex/NB01:/data zibra/zibra

### Build volume

Create a named data volume that mirrors local `build/` to `build/` within container:

    docker create --name zibra-build -v /Volumes/Meristem/build:/build zibra/zibra

This is to get data out of the Docker container. Note that the path to local directory has to be an absolute path.

### Start

Enter docker image:

    docker run -t -i --volumes-from zibra-data --volumes-from zibra-build zibra/zibra /bin/bash

Run single sample script within image:

    ./scripts/go_single_sample_r94.sh refs/KJ776791.2.fasta NB03 metadata/v2_500.amplicons.ver2.bed

## Results

  * [Initial coverage results from first library are here](depth-coverage/)
