# Zika bioinformatics pipeline -- Rhino

In this directory are all the scripts used in various parts of the bioinformatic analysis of USVI Zika analysis on the [Rhino cluster](https://github.com/blab/wiki/wiki/Rhino-cluster). For explanation of the Docker pipeline, see `docker_README.md`.

## Required modules
Different modules are required for different parts of the pipeline; they can be loaded using `module load <module name>` from Rhino.
- Barcoding:
  - `Python/2.7.13-foss-2016b-fh2`: Poretools
  - `Python/3.5.2-foss-2016b-fh1`: Porechop
- Basecalls:
  - __Albacore is not working at this time; issue with Rhino installation, it seems.__
- Pipeline:
  - `module load Python/3.6.1-foss-2016b-fh1`
  - `R/3.4.0-foss-2016b-fh1`
  - `BWA/0.7.15-foss-2016b`
  - `SAMtools/1.3.1-foss-2016b`
  - `nanopolish`


### `barcodes/`
Contains scripts and
