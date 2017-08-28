import time
import subprocess

BARCODES = [ 'NB%02d' % (s) for s in range(1,13) ]

def _get_prefix(wildcards):
    """Return the name that will be prepended onto final output consensus genomes
    as defined in the configuration file. If no prefix is given, a string of the
    current date.
    """
    if "prefix" in config:
        prefix = config["prefix"]
    else:
        prefix = time.strftime("%d-%m-%Y")

    return prefix

rule all:
    params:
        prefix=_get_prefix
    input:
        "{params.prefix}_good.fasta",
        "{params.prefix}_partial.fasta",
        "{params.prefix}_poor.fasta",
        ""

def _get_one_file():
    call = "find %s -name \"*.fast5\" | head -n 1" % config["raw_reads"]
    fname = subprocess.check_output(call,shell=True)
    return fname

ONE_FILE=_get_one_file()

def _get_raw_reads_dir(wildcards):
    return config["raw_reads"]

def _get_basecalled_reads_dir(wildcards):
    return config["basecalled_reads"]

def _get_albacore_config(wildcards):
    return config["albacore_config"]

rule basecall:
    params:
        one_file=ONE_FILE,
        raw_reads=_get_raw_reads_dir,
        basecalled_reads=_get_basecalled_reads_dir,
        cfg=_get_albacore_config
    input:
        raw="{params.raw_reads}/{params.one_file}"
    output:
        "{params.basecalled_reads}/workspace/{params.one_file}"
    shell:
        "read_fast5_basecaller.py -i {params.raw_reads} -t 8 --config {params.cfg} -r -s {params.basecalled_reads} -o fast5 -n 0"

def _get_demux_dir(wildcards):
    return config["demux_dir"]

rule extract_fasta:
    params:
        basecalled_reads=_get_basecalled_reads_dir,
        one_file=ONE_FILE
        demux_dir=_get_demux_dir
    input:
        "{params.basecalled_reads}/workspace/{params.one_file}"
    output:
        "{params.demux_dir}/nanopolish_full.fasta"
    shell:
        "nanopolish extract -b albacore -t template -o {output} {params.basecalled_reads}/workspace/"

rule demultiplex_full_fasta:
    params:
        demux_dir=_get_demux_dir
    input:
        "{params.demux_dir}/nanopolish_full.fasta"
    output:
        expand("{params.demux_dir}/{barcodes}.fasta.gz",barcodes=BARCODES)
    shell:
        "porechop -i {input} -b {params.demux_dir} --barcode_threshold 75 --threads 16 --check_reads 100000"

rule gunzip_demuxed_fastas:
    params:
        demux_dir=_get_demux_dir
    input:
        expand("{params.demux_dir}/{barcodes}.fasta.gz",barcodes=BARCODES)
    output:
        expand"({params.demux_dir}/{barcodes}.fasta",barcodes=BARCODES)
    shell:
        "gunzip -c {input} > {output}"

def _get_samples(wildcards):
    s = config['samples']
    samples = " ".join(s)
    return samples

rule pipeline:
    params:
        dimension=config['dimension']
        samples=_get_samples
    shell:
        "python pipeline/scripts/pipeline.py --samples {params.samples} --dimension {params.dimension}"
