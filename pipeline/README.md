# Zika bioinformatics pipeline

This directory contains scripts and files to run the bioinformatic analysis of Zika genomes from the USVI. Broadly, there are 3 distinct pipelines.

  1) Pipeline for genomes sequenced on the Illumina MiSeq.
  2) Pipeline for genomes sequenced on the MinION, that runs on the Rhino cluster at Fred Hutch.
  3) Pipeline for genomes sequenced on the MinION, run via a Docker image (currently deprecated).

The vast majority of this README describes the bioinformatic pipeline to analyze MinION data on the [Rhino cluster](https://github.com/blab/wiki/wiki/Rhino-cluster). Information on the MiSeq pipeline can be found in the [miseq_pipeline directory.](miseq_pipeline/)

For explanation of the Docker pipeline, see `docker_README.md`.

----------
## Running the pipeline from Nano (Bedford Lab personal Linux box):

_This pipeline uses the Conda package manager and Snakemake to compile the entire pipeline quickly on a laptop computer. One run of the pipeline will correspond to one sequencing library. Note that this changes may be necessary if run on an operating system other than Ubuntu 16.04 LTS._

#### Setting up and running the pipeline

1. Download Albacore

Download the `.deb` for Albacore from the [ONT Community](https://community.nanoporetech.com/downloads). If a version other than 2.0.2 is downloaded, modify the `dpkg` command in `install.sh` appropriately. Note that you need to have login credentials with Nanopore to access the Albacore software.

2. Download and install the pipeline from the github repo:
  ```
  git clone https://github.com/blab/zika-seq.git
  cd zika-seq
  git checkout create_wrap
```
If you are running Ubuntu 16.04:
```
  ./install_Ubuntu16_04.sh
```
If you are running Mac OSX:
```
  ./install_MacOSX.sh
```
  If `./<install_script>` fails, you may need to run `chmod 700 <install_script>` before rerunning.

Note: If you are running the install on Mac OSX, there's a possibility that nanopolish will not install properly. If so, you can do a manual install of nanopolish with:

`git clone --recursive https://github.com/jts/nanopolish.git`

`cd nanopolish`

`make CXX=g++-7 CC=gcc-7`

Then, add nanopolish to your $PATH with `export PATH=$PATH:~/your/path/to/nanopolish/`

We've found that the install tends to fail if you make without referencing these exact versions of g++ and gcc. If this still doesn't work, it might be an issue with Xcode install of CommandLineTools (we've found that this happens on OSX High Sierra). In this case, try running `xcode-select --install` to properly install CommandLineTools, then re-cloning and making the nanopolish directory.


3. Make any directories that you have specified in the config file that do not already exist. We usually put these directories in the repo itself. If it is your very first time running the repo you'll need to make `build` and `data` directories. For subsequent runs, you'll need to make library-specific directories that contain `basecalled_reads`, `process`, and `demux` directories.

4. Open `cfg.py` and change config information as appropriate:
  - `raw_reads` : directory containing un-basecalled `.fast5` numbered directories.
  - `dimension` : sequencing dimension (1d or 2d)
  - `demux_dir` : path to directory where demultiplexing will take place
  - `build_dir` : path to output location (`zika-seq/build`)
  - `samples` : list of all samples that are included for the library that will be processed
  - `albacore_config` : name of the config file to be used during basecalling by Albacore
  - `prefix` : prefix to prepend onto output consensus genome filenames

  Reminder: Make sure that all paths to directory paths listed in `demux_dir`, `build_dir` exist prior to running `snakemake`.

5. Run the pipeline:
  ```
  source activate zika-seq
  snakemake --use-conda
  ```

#### Tips

If this is your first time using `conda` and `snakemake`, you'll need to add the path to your version of miniconda to your bash profile. To check whether a path already exists, do `echo $PATH`. If there isn't a path to conda in the bash profile already, make one using `export PATH=$PATH:your/path/to/miniconda3/bin`.

Remember to move all final files that you want to keep to a separate directory, as everything in the `build` directory will be overwritten next time you run the pipeline (with a different config file).

If you had a single sample that failed, rather than re-running the pipeline on all the samples again, change the config file so that the pipeline will only run on the _single_ failed sample.

#### Inside the pipeline...

When Albacore basecalls the raw `fast5` reads, it makes a `workspace` directory, which then contains a `pass` and a `fail` directory. Both the `pass` and the `fail` directory contain a single `fastq`. The fastq in the `pass` directory contains the high quality reads, and the pipeline only uses this `fastq`.

The basecalled `fastq` file serves as input to [`porechop`](https://github.com/rrwick/Porechop), the program which performs barcode demultiplexing. Porechop writes a single `fastq` file for each barcode to the `demux` directory you created.

Next we run `pipeline.py`, which is a large script that references other custom python scripts and shell scripts to do run every other step of the pipeline. This is what is occurring in `pipeline.py`.

1.  Map sample IDs and run information that is contained in the `samples.tsv` and the `runs.tsv` files. This allows us to know which sample IDs are linked to which barcodes and primers etc.

2. Convert demuxed `fastq` files to `fasta` files. (If you look in the `demux` directory after a pipeline run you'll see both file types present for each barcode).

3. Using the linked sample and run information, combine the two barcode `fastq` or `fasta` files that represent pool 1 and pool 2 of the same sample. The pool-combined files are written as `<sample ID>.fasta` and `<sample ID>.fastq` to the `build` directory you made previously.

4. Map the reads in `<sample ID>.fasta` to a reference sequence using `bwa-mem`, and pipe to `samtools` to get a sorted bam file. Output files are labeled `<sample ID>.sorted.bam`.

5. Trim primer sequences out of the bam file using the custom script `align_trim.py`. Output files are labeled `<sample ID>.trimmed.sorted.bam`.

6. Use [`nanopolish`](https://github.com/jts/nanopolish) to call SNPs more accurately. This is a two step process. First step is using `nanopolish index` to create a map of how basecalled reads in `<sample ID>.fastq` are linked to the raw reads (signal level data). This step takes a while since for every sample all the program needs to iterate through all the raw reads.  Next, `nanopolish variants` will walk through the indexed `fastq` and a reference file, determine which SNPs are real given the signal level data, and write true SNPs to a VCF file.

7. Run `margin_cons.py` to walk through the reference sequence, the trimmed bam file, and the VCF file. This script looks at read coverage at a site, masking the sequence with 'N' if the read depth is below a hardcoded threshold (we use a threshold of 20 reads). If a site has sufficient coverage to call the base, either the reference base or the variant base (as recorded in the VCF) is written to the consensus sequence. Consensus sequences are written to `<sample ID>_complete.fasta`. The proportion of the genome that has over 20x coverage and over 40x coverage is logged to `<sample_ID>-log.txt`.

8. Consensus sequences are bundled together into `good`, `partial`, or `poor` files depending on the percent of the genome that was sequenced. Good quality genomes are >80% complete, partial genomes are between 50% and 80% complete, and poor genomes are <50% complete.


## Running the MinION bioinformatic pipeline on Rhino (Fred Hutch cluster)

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
  1. Load Poretools with `ml nanopolish/0.7.1-foss-2016b` (Uses Python 2).
  2. Change working directory to `<BASECALLED_READS_DIRECTORY>/workspace/demux/`.
  3. Submit job as `sbatch --time=96:00:00 --mem=32000 --mail-type=END,FAIL --mail-user=<EMAIL_ADDRESS> --wrap="$EBROOTNANOPOLISH/nanopolish extract -b albacore -t template -o nanopolish_full.fasta ../"`.
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

## Description of directories within `pipeline`:

#### `demux/`
Contains scripts and files relating to the demultiplexing of basecalled reads.

#### `basecall/`
Contains scripts used to prepare libraries for basecalling. __This directory is a work in progress, it will be made more user friendly soon.__

#### `scripts/`
Contains the majority of the scripts necessary to make consensus genomes from basecalled, demultiplexed reads, as well as for the generation of summary statistics and figures.
