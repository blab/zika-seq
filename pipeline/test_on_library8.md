# Library 8 full pipeline test -- 21 June 2017

## Data directories
Input: `/fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/raw_reads/`
Basecalls: `/fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/test/basecalled_reads/`

## Protocol

1. Submitted albacore:

```
bpotter@rhino1:/fh/fast/bedford_t/zika-seq$ module load Python/3.5.2-foss-2016b-fh2
bpotter@rhino1:/fh/fast/bedford_t/zika-seq$ sbatch --time=96:00:00 --mem=30000 --mail-type=END,FAIL --mail-user=bpotter@fhcrc.org --wrap="read_fast5_basecaller.py -i /fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/raw_reads/ -t 8 -c r94_450bps_linear.cfg -r -s /fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/test/basecalled_reads/ -o fast5"
Submitted batch job 52101562
```
Just over 2 days of Albacore runtime:

SLURM Job_id=52101562 Name=wrap Ended, Run time 2-00:50:09, COMPLETED, ExitCode 0

2. Moved all basecalled files into `workspace` directory:

```
bpotter@rhino1:/fh/fast/bedford_t/zika-seq$ python pipeline/basecall/prep_lib.py --library 8_test
```
3. Used Poretools to construct fasta of all basecalled reads:

```
bpotter@rhino3:/fh/fast/bedford_t/zika-seq$ ml Python/2.7.13-foss-2016b-fh2
bpotter@rhino3:/fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/test/basecalled_reads/workspace/demux$ sbatch --time=48:00:00 --mem=20000 --mail-type=END,FAIL --mail-user=bpotter@fhcrc.org --wrap="poretools fasta ../ > lib8-test-full.fasta"
Submitted batch job 52551608
```
