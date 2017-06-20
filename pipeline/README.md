# Zika bioinformatics pipeline -- Rhino

In this directory are all the scripts used in various parts of the bioinformatic analysis of USVI Zika analysis on the [Rhino cluster](https://github.com/blab/wiki/wiki/Rhino-cluster). For explanation of the Docker pipeline, see `docker_README.md`.

## Required modules
Different modules are required for different parts of the pipeline; they can be loaded using `module load <module name>` from Rhino.
- Barcoding:
  - `Python/2.7.13-foss-2016b-fh2`: Poretools
  - `Python/3.5.2-foss-2016b-fh1`: Porechop
- Basecalls:
  - `Python/3.5.2-foss-2016b-fh2`
  - __Albacore is not working at this time; issue with Rhino installation, it seems.__
- Pipeline:
  - `module load Python/3.6.1-foss-2016b-fh1`
  - `R/3.4.0-foss-2016b-fh1`
  - `BWA/0.7.15-foss-2016b`
  - `SAMtools/1.3.1-foss-2016b`
  - `nanopolish`


### `barcodes/`
Contains scripts and files relating to the demultiplexing of basecalled reads.

### `basecall/`
Contains scripts used to prepare libraries for basecalling. __This directory is a work in progress, it will be made more user friendly soon.__

### `scripts/`
Contains the majority of the scripts necessary to make consensus genomes from basecalled, demultiplexed reads, as well as for the generation of summary statistics and figures.
