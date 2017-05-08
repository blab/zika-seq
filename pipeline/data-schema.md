# Data schema

## Nanopore reads

Input data to the Zika pipeline arrives in the `data/` directory.

  - `data/`
    - `usvi-library1-2016-12-10/` - library
      - `raw_reads/` - squiggle graphs in fast5 format; all subdirectores are written automatically by MinKNOW
        - `0/` - contains ~4000 raw `.fast5` files
        - `1/` - contains ~4000 raw `.fast5` files
        - etc ...
      - `basecalled_reads/` - basecalled with Albacore 1.0.4; all subdirectores are written automatically by Albacore
        - `sequencing_summary.txt` - summary of Albacore run; automatically made by Albacore
        - `pipeline.log` - summary of Albacore run; automatically made by Albacore
        - `workspace` - contains all basecalled, demultiplexed reads
          - `barcode01/` - basecalled, demultiplexed reads with ONT barcode NB01
            - `0/` - contains ~4000 basecalled, demultiplexed `.fast5` files
            - `1/` - contains ~4000 basecalled, demultiplexed `.fast5` files
            - etc ...
          - `barcode02/` - basecalled, demultiplexed reads with ONT barcode NB02
            - `0/` - contains ~4000 basecalled, demultiplexed `.fast5` files
            - `1/` - contains ~4000 basecalled, demultiplexed `.fast5` files
            - etc ...
          - ...
          - `barcode12/` - basecalled, demultiplexed reads with ONT barcode NB12
            - `0/` - contains ~4000 basecalled, demultiplexed `.fast5` files
            - `1/` - contains ~4000 basecalled, demultiplexed `.fast5` files
            - etc ...
          - `unclassified/` - basecalled, non-demultiplexed reads
            - `0/` - contains ~4000 basecalled, non-demultiplexed `.fast5` files
            - `1/` - contains ~4000 basecalled, non-demultiplexed `.fast5` files
            - etc ...

## Sample metadata

Sample metadata for the Zika pipeline arrives in the `samples/` directory.

  - `samples/`
    - `samples.tsv` - line list of sample metadata
    - `runs.tsv` - line list of run metadata

### `samples.tsv`

Must be `tsv` formatted. Keyed off of column headers rather than column order.

sample_id | strain  | collection_date | country | division     | location
--------- | ------- | --------------- | ------- | ---------- | ---------
ZBRD116   | ZBRD116 | 2015-08-28      | brazil  | alagoas    | arapiraca
ZBRC301   | ZBRC301 | 2015-05-13      | brazil  | pernambuco | paulista

### `runs.tsv`

Must be `tsv` formatted. Keyed off of column headers rather than column order.

run_name | barcode_id | sample_id | primer_scheme
-------- | ---------- | --------- | -------------
library1 | NB01       | ZBRD116   | v2_500.amplicons.ver2
