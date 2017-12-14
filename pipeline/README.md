# Zika bioinformatics pipeline

This directory contains scripts and files to run the bioinformatic analysis of Zika genomes from the USVI. Broadly, there are 3 distinct pipelines.

  1) Pipeline for genomes sequenced on the Illumina MiSeq.
  2) Pipeline for genomes sequenced on the MinION, that runs on the Rhino cluster at Fred Hutch.
  3) Pipeline for genomes sequenced on the MinION, run via a Docker image (currently deprecated).

The vast majority of this README describes the bioinformatic pipeline to analyze MinION data on the [Rhino cluster](https://github.com/blab/wiki/wiki/Rhino-cluster). Information on the MiSeq pipeline can be found in the [miseq_pipeline directory.](miseq_pipeline/)

For explanation of the Docker pipeline, see `docker_README.md`.

----------

## Running the MinION bioinformatic pipeline on Rhino

### Required modules
Different modules are required for different parts of the pipeline; they can be loaded using `module load <module name>` from Rhino. Descriptions of when each module should be loaded are in the step by step instructions below.

- Basecalling:
  - `Python/3.5.2-foss-2016b-fh2`: Albacore 1.2.1
- Barcode demultiplexing:
  - `Python/2.7.13-foss-2016b-fh2`: Poretools
  - `Python/3.5.2-foss-2016b-fh1`: Porechop
- Consensus genome generation: (all necessary for `pipeline.py` step)
  - `module load Python/3.6.1-foss-2016b-fh1`
  - `R/3.4.0-foss-2016b-fh1`
  - `BWA/0.7.15-foss-2016b`
  - `SAMtools/1.3.1-foss-2016b`
  - `nanopolish/0.7.1-foss-2016b`: __Breaks one step of pipeline, waiting for update from SciComp.__

### Step-by-step running instructions

##### Notes
* <RAW_READS_DIRECTORY> needs to be a global path, for instance `/fh/fast/bedford_t/zika-seq/data/<RUN>/raw_reads/`
* <BASECALLED_READS_DIRECTORY> is assumed to be `/fh/fast/bedford_t/zika-seq/data/<RUN>/alba121/` by `pipeline.py`. If this is not the case, errors will arise on step 4.2.
* <CONFIG_FILE> is `r94_250bps_nsk007_2d.cfg` for 2D reads and `r94_450bps_linear.cfg` for 1D reads.
* <SAMPLES_TO_RUN> are, for example: `VI1 VI2 VI3`.

##### Basecall raw reads using Albacore 1.2.1 installed on the cluster:
  1. Load Albacore with `module load Python/3.5.2-foss-2016b-fh2`.
  2. Sbatch out the Albacore job as follows `sbatch --time=48:00:00 --mem=20000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="read_fast5_basecaller.py -i <RAW_READS_DIRECTORY> -t 8 --config <CONFIG_FILE> -r -s <BASECALLED_READS_DIRECTORY> -o fast5`

##### Reorganizing the `workspace` directory:

Albacore writes basecalled `fast5` files into subdirectories within the `workspace` directory. Each directory is given a number (starts at `0/`). Albacore writes 4000 files into a directory, and then makes a new directory, in numerical order. In the end, you'll have a lot of subdirectories within workspace. Importantly, [Porechop](https://github.com/rrwick/Porechop), which we use for demultiplexing the barcoded reads, requires that all the files be in a single directory. Therefore, the first step once you have basecalled reads is to use [`prep_lib.py`](basecall/prep_lib.py) to move all the files out of the subdirectories into a single directory. This script also makes a directory called `demux/`, which is where Porechop will write the demultiplexed fasta files.

   1. Submit the library prep job as `sbatch --time=48:00:00 --mem=10000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="python pipeline/basecall/prep_lib.py --library <LIBRARY_NUMBER>`.

###### Note: Your library might be huge, in which case there will be too many files in `workspace` for Porechop to handle.

If this is the case, then you will need to split the files back up into smaller subdirectories, and you'll need to submit a Porechop job _for each subdirectory_. If your library has more than 4 million reads, you can use the [`split_1d_library.py`](demux/split_1d_library.py) to sort reads into directories that contain 500,000 reads per directory.

##### Extract fasta files from basecalled reads and demultiplex reads based on barcoding:
  1. Load Poretools with `module load Python/2.7.13-foss-2016b-fh2` (Uses Python 2).
  2. Change working directory to `<BASECALLED_READS_DIRECTORY>/workspace/demux/`.
  3. Submit job as `sbatch --time=48:00:00 --mem=20000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="poretools fasta ../ > <FILENAME.fasta>`.
  4. Load Porechop with `module load Python/3.5.2-foss-2016b-fh1` (uses Python 3).
  5. Submit job as `sbatch --time=24:00:00 --mem=20000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="porechop -i <FILENAME.fasta> -b . --barcode_threshold 75 --threads 16 --check_reads 100000"`.
  6. Once job completes, run `gunzip NB*` from within the directory to unzip the files in preparation for consensus genome generation.

Porechop will write out a fasta file for each barcode in the `demux/` directory (for example, `demux/NB01.fasta`). If you had a large library, and you split basecalled reads into subdirectories of 500,000 files each, and submitted a Porechop job for _each_ subdirectory, then you'll have multiple fastas with the demultiplexed reads. In such cases, use [`cat_demux_fastas.py`](demux/cat_demux_fastas.py) to consolidate all of the reads into a single fasta for each barcode.

##### Generating consensus sequences:

   1. Change working directory to `/fh/fast/bedford_t/zika-seq/`
   2. Load all of the following modules (`ml` is a default Rhino alias for `module load`):
   - `ml R/3.4.0-foss-2016b-fh1`   
   - `ml Python/3.6.1-foss-2016b-fh1`
   - `ml BWA/0.7.15-foss-2016b`
   - `ml SAMtools/1.3.1-foss-2016b`
   - `ml nanopolish/0.7.1-foss-2016b`
   3. Submit the job as `sbatch --time=48:00:00 --mem=30000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="python pipeline/scripts/pipeline.py --samples <SAMPLES_TO_RUN> --dimension <DIMENSION>"`.

The [`pipeline.py`](scripts/pipeline.py) script does all the heavy lifting in terms of alignment to a reference, calling variants, and writing consensus genomes. Details are in a separate [`README`](scripts/README.md). Note that output is written to the `build/` directory.

## Running the pipeline from Nano:
_This pipeline uses the Conda package manager and Snakemake to compile the entire pipeline quickly on a laptop computer. One run of the pipeline will correspond to one sequencing library. Note that this changes may be necessary if run on an operating system other than Ubuntu 16.04 LTS._

1. Download Albacore

Download the `.deb` for Alabacore from the (ONT)[https://nanoporetech.com/]. If a version other than 2.0.2 is downloaded, modify the `dpkg` command in `install.sh` appropriately.

2. Download and install from the github repo:
```
git clone https://github.com/blab/zika-seq.git
git checkout create_wrap
chmod 700 install.sh
./install.sh
```

3. Open `cfg.py` and change config information as appropriate:
  - `raw_reads` : directory containing un-basecalled `.fast5` numbered directories.
  - `dimension` : sequencing dimension (1d or 2d)
  - `demux_dir` : path to directory where demultiplexing will take place
  - `build_dir` : path to output location (`zika-seq/build`)
  - `samples` : list of all samples that are included for the library that will be processed
  - `albacore_config` : name of the config file to be used during basecalling by Albacore
  - `prefix` : prefix to prepend onto output consensus genome filenames

  Important: Make sure that all paths to directory paths listed in `demux_dir`, `build_dir` exist prior to running `snakemake`.

4. Run the pipeline:
  ```
  source activate zika-seq
  snakemake --use-conda
  ```

## Description of directories within `pipeline`:

#### `demux/`
Contains scripts and files relating to the demultiplexing of basecalled reads.

#### `basecall/`
Contains scripts used to prepare libraries for basecalling. __This directory is a work in progress, it will be made more user friendly soon.__

#### `scripts/`
Contains the majority of the scripts necessary to make consensus genomes from basecalled, demultiplexed reads, as well as for the generation of summary statistics and figures.
