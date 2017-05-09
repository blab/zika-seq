# Zika-USVI Pipeline
Pipeline for processing raw nanopore-sequenced Zika samples to create consensus genomes. Developed for use on Rhino.

## Requirements
Prior to running different parts of the pipeline, make sure that the requisite modules have been loaded on Rhino using `module load <module>`.

Confirm that the correct modules have been loaded using `module list`.

#### Basecalling:
* ONT Albacore 1.0.4; loaded with `module load Python/3.5.2-foss-2016b-fh1`

#### Consensus genomes:
* Python: `module load Python/3.6.1-foss-2016b-fh1`
* R: `module load R/3.4.0-foss-2016b-fh1`
* Nanopolish: `module load nanopolish`
  - Note: nanopolish must be run with command `$EBROOTNANOPOLISH/nanopolish`
* BWA: `module load BWA/0.7.15-foss-2016b`
* SAMtools: `module load SAMtools/1.3.1-foss-2016b`

## Running the pipeline

#### Basecalling:

Basecalling can be distributed to Rhino automatically using the command:

`python distribute_albacore.py -i </path/to/library/raw_reads> -o </path/to/library/basecalled_reads/> --email <user email>`

This command will distribute an SBATCH command for each of the subdirectores of `raw_reads/`, each containing roughly 4000 `.fast5` files. Emails will be sent to the given address for each subdirectory that failed, typically ~5-10 per distribution. Most failures happen at the beginning of the runs.

`distribute_albacore.py` has the following arguments that can be specified by the user (`inpath`, `outpath`, and `email` are required):

* `-i`, `--inpath`: path to input directory containing non-basecalled reads
  - Note this directory name _should not_ be followed by a terminal `/`
* `-o`, `--outpath`: path to output directory for basecalled reads
  - Note this directory name _should_ be followed by a terminal `/`
* `--dimension`: dimension of sequenced library; `1d` or `2d`
  - Default is `1d`
* `--email`: email address for sbatch notifications
* `--test`: make sure that albacore is working correctly
  - Passes albacore help call, `read_fast5_basecaller.py -h`, to SLURM. Should give a `slurm.out` file within a few minutes that contains the help messages for albacore. If there is an error with albacore (usually incorrect modules have been loaded, see above), it will show here.
* `--dryrun`: print commands to be run, but do not call them
  - Will print all SBATCH calls to the screen, but will not execute them. Usefull for making sure that all input is correct before submitting hundreds of jobs to SLURM.

It generally works run in 3 steps, to avoid hundreds or thousands of error emails being sent immediately:
1. Full command with `--test` to make sure that the proper modules have been loaded.
2. Full command with `--dryrun` to make sure that all input is correct. A common error happens if `--inpath` is specified with a terminal `/` (i.e. `path/to/library/raw_reads/` rather than `path/to/library/raw_reads`).
3. Full command without any test flags, to submit jobs to SLURM.

If multiple libraries are being submitted at the same time, it is recommended to allow 10-15 minutes to pass between `distribute_albacore` calls. Most errors happen during the first few minutes after submission, so separating calls by some time helps when determining which error email corresponded to which library.
