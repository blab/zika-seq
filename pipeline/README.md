# Zika bioinformatics pipeline -- Rhino

In this directory are all the scripts used in various parts of the bioinformatic analysis of USVI Zika analysis on the [Rhino cluster](https://github.com/blab/wiki/wiki/Rhino-cluster). For explanation of the Docker pipeline, see `docker_README.md`.

## Required modules
Different modules are required for different parts of the pipeline; they can be loaded using `module load <module name>` from Rhino.

- Basecalls:
  - `Python/3.5.2-foss-2016b-fh2`: Albacore 1.2.1
- Barcoding:
  - `Python/2.7.13-foss-2016b-fh2`: Poretools
  - `Python/3.5.2-foss-2016b-fh1`: Porechop
- Pipeline: (all necessary for `pipeline.py` step)
  - `module load Python/3.6.1-foss-2016b-fh1`
  - `R/3.4.0-foss-2016b-fh1`
  - `BWA/0.7.15-foss-2016b`
  - `SAMtools/1.3.1-foss-2016b`
  - `nanopolish`: __Breaks one step of pipeline, waiting for update fro SciComp.__

## To run
* <RAW_READS_DIRECTORY> is generally `/fh/fast/bedford_t/zika-seq/data/<RUN>/raw_reads/`
* <BASECALLED_READS_DIRECTORY> is assumed to be `/fh/fast/bedford_t/zika-seq/data/<RUN>/alba121/` by `pipeline.py`. If this is not the case, errors will arise on step 4.2.
* <CONFIG_FILE> is `r94_250bps_nsk007_2d.cfg` for 2D reads and `r94_450bps_linear.cfg` for 1D reads.
* <SAMPLES_TO_RUN> are, for example: `VI1 VI2 VI3`.
1. Barcoding:
   1. `sbatch --time=48:00:00 --mem=20000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="read_fast5_basecaller.py -i <RAW_READS_DIRECTORY> -t 8 --config <CONFIG_FILE> -r -s <BASECALLED_READS_DIRECTORY> -o fast5`
2. Workspace directory preparation:
   1. `sbatch --time=48:00:00 --mem=10000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="python pipeline/basecall/prep_lib.py --library <LIBRARY_NUMBER>`
3. FASTA extraction and barcoding:
   1. Change working directory to `<BASECALLED_READS_DIRECTORY>/workspace/demux/`
   2. `sbatch --time=48:00:00 --mem=20000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="poretools fasta ../ > <FILENAME.fasta>`
   3. `sbatch --time=24:00:00 --mem=20000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="porechop -i <FILENAME.fasta> -b . --barcode_threshold 75 --threads 16 --check_reads 100000"`
4. Pipeline:
   1. Change working directory to `/fh/fast/bedford_t/zika-seq/`
   2. `sbatch --time=48:00:00 --mem=20000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="python pipeline/scripts/pipeline.py --samples <SAMPLES_TO_RUN> --dimension <DIMENSION>`

### `barcodes/`
Contains scripts and files relating to the demultiplexing of basecalled reads.

### `basecall/`
Contains scripts used to prepare libraries for basecalling. __This directory is a work in progress, it will be made more user friendly soon.__

### `scripts/`
Contains the majority of the scripts necessary to make consensus genomes from basecalled, demultiplexed reads, as well as for the generation of summary statistics and figures.
