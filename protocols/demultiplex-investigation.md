#Overview
We want to find out what issues may be causing inconsistency in the basecalling and demultiplexing of our 1d and 2d reads, and we seek the best pipeline for getting the maximum number of good quality reads from both 1d and 2d libraries. In particular, we believe that Albacore results in lower quality basecalling on 2d libraries than on 1d libraries, which seems to cause large differences when we try to demultiplex.

During this, we find that there are huge differences in the quality of `fasta` generation when using `nanopolish extract` vs `poretools fasta`.

## Metrichor vs Albacore basecalling and demux
To test the differences in basecall/demux, we wanted to test the quality of the native Albacore demultiplexer compared to Porechop (recommended to us by ONT and Zibra) and the Metrichor demuxer. We also want to look at the difference in basecalling between Metrichor CLI (2d reads) and Albacore (1d and 2d reads).

| | Basecaller | Demuxer | Library (dimension) |
|-|:-----------|:--------|:--------------------|
|1| Metrichor | Metrichor | USVI Library 1 (2D) |
|2| Metrichor | Porechop | USVI Library 1 (2D) |
|3| Albacore | Porechop | USVI Library 1 (2D) |
|4| Albacore | Albacore | USVI Library 1 (2D) |
|5| Albacore | Albacore | USVI Library 7 (1D) |
|6| Albacore | Porechop | USVI Library 7 (1D) |

Metrichor basecalled (and demultiplexed) fast5 files already existed from initial processing of Library 1.
Albacore basecalls:
- `read_fast5_basecaller.py -i <INDIR> -t 8 -c <CONFIG> -r -s <OUTDIR> -o fast5 (--barcoding)`
  - This writes from the `-i` directory to the `-s` directory in `OUTDIR/workspace/barcodeXX/` if `--barcoding` flag is used and `OUTDIR/workspace/` if not
  - CONFIG is `r94_450bps_linear.cfg` for 1d libraries and `r94_250bps_2d.cfg` for 2d
  - Library 1:
    - INDIR: `/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/raw_reads/`
    - OUTDIR (barcoding): `/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/basecalled_reads/` _Note, this directory has been un-demultiplexed since running, so barcodeXX directories no longer exist_
    - OUTDIR (no barcoding): `/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/test/`
  - Library 7:
    - INDIR: `/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/raw_reads/`
    - OUTDIR (barcoding): `/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/basecalled_reads/` _Note, this directory has been un-demultiplexed since running, so barcodeXX directories no longer exist_
    - OUTDIR (no barcoding): `/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/test/`

### Library 1 (2D) -- Metrichor basecall
#### Metrichor demux
```
barcode01: 24415
barcode02: 21219
barcode03: 8897
barcode04: 31576
barcode05: 33059
barcode06: 215
barcode07: 19349
barcode08: 16519
barcode09: 12041
barcode10: 18536
barcode11: 25150
barcode12: 100
fail: 419911
```
This is 211076 classified out of 630987 total or 33% classified.
#### Porechop demux
```
Barcode  Reads   Bases       File           
NB01     17,531   7,658,568  ./NB01.fasta.gz
NB02     16,992   7,170,800  ./NB02.fasta.gz
NB03      6,022   2,540,559  ./NB03.fasta.gz
NB04     20,506   8,835,581  ./NB04.fasta.gz
NB05     26,015  11,437,361  ./NB05.fasta.gz
NB07     12,409   5,459,126  ./NB07.fasta.gz
NB08     11,551   5,023,361  ./NB08.fasta.gz
NB09      8,308   3,679,595  ./NB09.fasta.gz
NB10     12,316   5,428,981  ./NB10.fasta.gz
NB11     16,781   7,338,385  ./NB11.fasta.gz
NB12          1          35  ./NB12.fasta.gz
none     91,818  40,682,531  ./none.fasta.gz
```
This is 148432 classified out of 630987 total or 24% classified.
### Library 1 (2D) -- Albacore basecall
#### Albacore demux
```
barcode01: 440
barcode02: 399
barcode03: 183
barcode04: 491
barcode05: 578
barcode06: 17
barcode07: 434
barcode08: 387
barcode09: 268
barcode10: 392
barcode11: 600
barcode12: 4
unclassified: 606755
```
This is 4193 classified out of 610948 total or 0.7% classified.
#### Porechop demux
```
Barcode  Reads   Bases      File           
NB01        102     44,556  ./NB01.fasta.gz
NB02         97     40,426  ./NB02.fasta.gz
NB03         39     16,593  ./NB03.fasta.gz
NB04         79     35,008  ./NB04.fasta.gz
NB05        134     59,450  ./NB05.fasta.gz
NB06          5        701  ./NB06.fasta.gz
NB07         78     35,250  ./NB07.fasta.gz
NB08        101     43,031  ./NB08.fasta.gz
NB09         69     30,320  ./NB09.fasta.gz
NB10         81     36,436  ./NB10.fasta.gz
NB11        128     56,433  ./NB11.fasta.gz
none     15,592  5,208,138  ./none.fasta.gz
```
This is 913 classified out of 610948 total or 0.1% classified.
### Library 7 (1D)
#### Albacore demux
```
barcode01: 325651
barcode02: 334437
barcode03: 344680
barcode04: 351069
barcode05: 422758
barcode06: 2934
barcode07: 281634
barcode08: 277140
barcode09: 358521
barcode10: 14
barcode11: 359125
barcode12: 360467
unclassified: 1762939
```
This is 3418430 classified out of 5181369 total or 66% classified.
#### Porechop demux
**_Still Running_**

## Nanopolish vs Poretools
We also wanted to test the quality of `fasta` extraction by [Nanopolish](https://github.com/jts/nanopolish) and [Poretools](https://github.com/arq5x/poretools). For this, we took 1000 `fast5` subsets of Libraries 1 and 7, and ran either `nanopolish extract` or `poretools fasta` to make four FASTAs, then counted the number of reads represented in the FASTA using `grep -c '>' <FASTA>`. Those read counts were:

| | Library 1 | Library 7 |
|-|-----------|-----------|
| Nanopolish | 26 | 997 |
| Poretools | 1,045 | 1,026 |

Both sets generated using Albacore -> fast5 files detailed above
* 1d: `/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/test/workspace/demux/poretools_test`
* 2d: `/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/test/demux/poretools_test`
