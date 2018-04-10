# Zika-USVI Pipeline
Pipeline for processing raw nanopore-sequenced Zika samples to create consensus genomes. Developed for use on Rhino.

## `pipeline.py`

This script constitutes the bulk of the pipeline, generating consensus genomes from demultiplexed FASTAs.

### Arguments:
- `--data_dir`: Directory containing all libraries and data. Default is `/fh/fast/bedford_t/zika-seq/data/`. Contained in this directory should be demultiplexed FASTAs in the subfolder `<run>/alba121/workspace/demux/`
- `--samples_dir`: Directory containing a `runs.tsv` file and a `samples.tsv` file. These files are parsed by `pipeline.py` to generate metadata for each sample. Default is `/fh/fast/bedford_t/zika-seq/samples/`.
- `--build_dir`: Directory into which all intermediary files and output from `pipeline.py` will be written. Default is `/fh/fast/bedford_t/zika-seq/build/`.
- `--prefix`: String that will be prepended onto all output consensus genome files. Default is `ZIKA_USVI`.
- `--samples`: Names of samples that should be processed. Acceptable samples can be found in `<samples_dir>/samples.tsv`. Samples should be listed separated by spaces. If excluded from command, all samples listed in `samples.tsv` will be processed.
- `--dimension`: Dimension of the library being processed. Options are `1d` or `2d`; default is `2d`.
- `--run_steps`: Numbered steps to run (explained below in more detail):
  1. Construct sample FASTAs
  2. Process sample FASTAs
  3. Gather consensus FASTAs
  4. Generate coverage overlap plots
  5. Calculate per-base error rates

### Pipeline overview
##### Construct sample FASTAs
FASTAs for each sample are constructed by concatenating the two demultiplexed FASTAs that correspond to a given sample. The complete FASTA is written to `<build_dir>`.
##### Process sample FASTAs
Take a complete FASTA and construct a consensus genome. This is done by calling the script `fasta_to_consensus_<dimension>.sh`, which does the following:
  1. Aligns reads with `bwa mem`.
  2. Trims alignment to primer start sites with `align_trim.py`.
  3. Normalizes reads to sequencing depth of 500 in order to save time with `align_trim.py`.
  4. Creates a sorted BAM file with `samtools sort` and `samtools index`.
  5. Variant calling using `nanopolish variants`.
  6. Consensus genome construction with `margin_cons.py`.
Output consensus genomes are written to `<build_dir>` by sample name.
##### Gather consensus FASTAs
Iterates over all samples in `<build_dir>` to determine the percent of the genome that was called as 'N', and joins the consensus FASTAs into one of three files depending on percent N:
  - `ZIKA_USVI_good.fasta`: Less than 20% N
  - `ZIKA_USVI_partial.fasta`: Between 20% and 50% N
  - `ZIKA_USVI_poor.fasta`: Less than 50% N
##### Generate coverage overlap plots
Looks at the BAM files to determine the depth of sequencing at each site along the reference genome. Using this, plots are generated using `depth_coverage.R`. PNG files for each plot are written to `<build_dir>`.
* _Note:_ These plots can be generated without consensus genomes by running `depth_process.py` before calling `depth_coverage.R`.
##### Calculate per-base error rates
Walks through VCF files made by `nanopolish variants` to determine per-base error rates. _Currently broken_
