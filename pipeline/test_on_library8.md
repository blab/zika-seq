# Library 8 full pipeline test -- 21 June 2017

## Data directories
Input: `/fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/raw_reads/`
Basecalls: `/fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/test/basecalled_reads/`

## Protocol

### Submitted albacore (12:10):
```
bpotter@rhino1:/fh/fast/bedford_t/zika-seq$ sbatch --time=96:00:00 --mem=30000 --mail-type=END,FAIL --mail-user=bpotter@fhcrc.org --wrap="read_fast5_basecaller.py -i /fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/raw_reads/ -t 8 -c r94_450bps_linear.cfg -r -s /fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/test/basecalled_reads/ -o fast5"
Submitted batch job 52101562
```
