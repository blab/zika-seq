import time
import subprocess
from cfg import config

BARCODES = [ 'NB%02d' % (s) for s in range(1,13) ]
DEMUX_DIR=config['demux_dir']
BASECALLED_READS=config['basecalled_reads']
RAW_READS = config["raw_reads"]

def get_one_file():
    call = "find %s -name \"*.fast5\" | head -n 1" % config["raw_reads"]
    fname = subprocess.check_output(call,shell=True)
    if type(fname) != str:
        fname = str(fname)[2:-3]
    fname = fname.replace(config["raw_reads"],"")[1:]
    print("Example file: %s" % (fname))
    return fname

ONE_FILE=get_one_file()

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

def _get_albacore_config(wildcards):
    return config["albacore_config"]

rule basecall:
    params:
        cfg=_get_albacore_config
    input:
        raw="%s/%s" % (RAW_READS, ONE_FILE)
    output:
        "%s/workspace/pipeline.log" % (BASECALLED_READS)
    shell:
        "read_fast5_basecaller.py -i %s -t 8 --config {params.cfg} -r -s %s -o fast5 -n 0" % (RAW_READS, BASECALLED_READS)

rule extract_fasta:
    input:
        "%s/workspace/pipeline.log" % (BASECALLED_READS)
    output:
        "%s/nanopolish_full.fasta" % (DEMUX_DIR)
    shell:
        "nanopolish extract -b albacore -t template -o {output} %s/workspace/" % (BASECALLED_READS)

rule demultiplex_full_fasta:
    input:
        "%s/nanopolish_full.fasta" % (DEMUX_DIR)
    output:
        "%s/{barcodes}.fasta.gz" % (DEMUX_DIR)
    shell:
        "porechop -i {input} -b %s --barcode_threshold 75 --threads 16 --check_reads 100000" % (DEMUX_DIR)

rule gunzip_demuxed_fastas:
    input:
        "%s/{barcodes}.fasta.gz" % (DEMUX_DIR)
    output:
        "%s/{barcodes}.fasta" % (DEMUX_DIR)
    shell:
        "gunzip -c {input} > {output}"

def _get_samples(wildcards):
    """ Build a string of all samples that will be processed in a pipeline.py run.
    """
    s = config['samples']
    samples = " ".join(s)
    return samples

rule pipeline:
    params:
        prefix=_get_prefix,
        dimension=config['dimension'],
        samples=_get_samples
    input:
        expand("%s/{barcodes}.fasta" % (DEMUX_DIR),barcodes=BARCODES)
    output:
        "{params.prefix}_good.fasta",
        "{params.prefix}_partial.fasta",
        "{params.prefix}_poor.fasta"
    shell:
        "python pipeline/scripts/pipeline.py --samples {params.samples} --dimension {params.dimension}"
