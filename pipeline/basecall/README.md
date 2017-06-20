# Basecalls
This directory contains scripts related to preparing raw `.fast5` reads for basecalling.
This involves moving files to consistent locations, normalizing read counts for uniform computational times, and generating summary statistics on all or a subsets of libraries.

_Note that this directory does not handle any reads from USVI library 2, as that library was sequenced using miSeq, not minION._

### Sbatch job submission
For all operations involving more than 10,000 files (everything other than `grab_1000_fast5s.py`), it is recommended to use Slurm to hand the job off to Rhino, rather than running in a terminal window. Our default Slurm command is:

```
sbatch --time=24:00:00 --mem=10000 --mail-type=END,FAIL --mail-user=<username>@fhcrc.org --wrap="<job call>"
```

This specifies 24h of compute time, and 10GB of memory for the job. Sometimes, increasing `--mem` or `--time` is necessary. This also sends an email when the job either fails or ends successfully. The `--wrap` command takes the exact command that would normally be run from the terminal, and it inherits whichever modules are loaded at the time that the job is submitted.

Many of the scripts in this directory can be used to hand off commands that move tons of files to Slurm, so that they are not reliant on keeping a terminal window open the whole time that they run.

### Scripts

- `clear_workspace.py`: deletes the `basecalled_reads/workspace` directory for the chosen library
- `count_ratios.py`: Looks at all libraries and quantifies the % successful demux by Albacore's native demuxer
- `distribute_albacore.py`: For 1d libraries, hands Albacore jobs to slurm for every subdirectory in the `raw_reads` folder
- `grab_1000_fast5s.py`: Moves 1000 fast5 files (either 1d or 2d) into a test directory for small-scale analysis
- `prep_lib.py`: Moves files out of numbered subdirectories of `basecalled_reads/workspace`, deletes those numbered folders, and creates a `demux` folder for downstream steps
- `remove_tmp.py`: Used for handling library8 reads that did not initially sync to Rhino correctly during minKNOW run
- `subdivide.py`: un-does `prep-lib.py`; used for library8 after all reads were successfully copied to hand multiple jobs to slurm
