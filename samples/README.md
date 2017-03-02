# Sample metadata

Sample metadata for the Zika pipeline arrives in the `samples/` directory. This should be mounted to `samples/` in the Docker container.

  - `samples/`
    - `samples.tsv` - line list of sample metadata
    - `runs.tsv` - line list of run metadata

### `samples.tsv`

Must be `tsv` formatted. Keyed off of column headers rather than column order. `sample_id` is keyed to `accession` in [fauna](https://github.com/nextstrain/fauna) because these lack Genbank accessions.

sample_id | strain      | collection_date | country | division       | location    | usvi_sample_id | seq_platform |
--------- | ----------- | --------------- | ------- | -------------- | ----------- | -------------- | ------------ |
VI1       | USVI/1/2016 |	2016-09-28      | usvi    | saint_croix    | saint_croix | 16VI2395U      | minion       |

### `runs.tsv`

Must be `tsv` formatted. Keyed off of column headers rather than column order.

run_name                 | barcode_id | sample_id | primer_scheme
------------------------ | ---------- | --------- | -------------
usvi-library1-2016-12-10 | NB01       | VI1       | v2_500.amplicons.ver2
