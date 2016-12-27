# Data organization

  - `libraries`
    - `usvi-library1-2016-12-10`
      - `raw_reads`: FAST5 reads before Metrichor processing
      - `basecalled_reads`: FAST5 reads after Metrichor processing
        - `pass_demultiplex`:
          - `NB01`, etc...: List of `.fast5` files
        - `fail`: List of `.fast5` files
