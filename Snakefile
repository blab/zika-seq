import time
import subprocess
from cfg import config

BARCODES = [ 'BC%02d' % (s) for s in range(1,13) ]
DEMUX_DIR=config['demux_dir']
BASECALLED_READS=config['basecalled_reads']
RAW_READS = config["raw_reads"]
BUILD_DIR = config["build_dir"]
PREFIX = config["prefix"]

def get_one_file():
    call = "find %s -name \"*.fast5\" | head -n 1" % config["raw_reads"]
    fname = subprocess.check_output(call,shell=True)
    if type(fname) != str:
        fname = str(fname)[2:-3]
    fname = fname.replace(config["raw_reads"],"")[1:]
    print("Example file: %s" % (fname))
    return fname

ONE_FILE=get_one_file()

rule all:
    params:
        prefix=PREFIX,
        build=BUILD_DIR
    input:
        "%s/%s_good.fasta" % (BUILD_DIR, PREFIX),
        "%s/%s_partial.fasta" % (BUILD_DIR, PREFIX),
        "%s/%s_poor.fasta" % (BUILD_DIR, PREFIX)

def _get_albacore_config(wildcards):
    return config["albacore_config"]

rule basecall:
    params:
        cfg=_get_albacore_config
    input:
        raw="%s/%s" % (RAW_READS, ONE_FILE)
    output:
        "%s/pipeline.log" % (BASECALLED_READS)
    shell:
        "read_fast5_basecaller.py -i %s -t 8 --config {params.cfg} -r -s %s -o fast5 -n 0" % (RAW_READS, BASECALLED_READS)

rule extract_fasta:
    input:
        "%s/pipeline.log" % (BASECALLED_READS)
    output:
        "%s/nanopolish_full.fasta" % (DEMUX_DIR)
    conda:
        "envs/anaconda.nanopolish-env.yaml"
    shell:
        "nanopolish extract -b albacore -t template -o {output} %s/workspace/pass" % (BASECALLED_READS)

rule demultiplex_full_fasta:
    input:
        "%s/nanopolish_full.fasta" % (DEMUX_DIR)
    output:
        "%s/{barcodes}.fasta" % (DEMUX_DIR)
    shell:
        "porechop -i {input} -b %s --barcode_threshold 75 --threads 16 --check_reads 100000" % (DEMUX_DIR)

# rule gunzip_demuxed_fastas:
#     input:
#         "%s/{barcodes}.fasta.gz" % (DEMUX_DIR)
#     output:
#         "%s/{barcodes}.fasta" % (DEMUX_DIR)
#     shell:
#         "gunzip -c {input} > {output}"

def _get_samples(wildcards):
    """ Build a string of all samples that will be processed in a pipeline.py run.
    """
    s = config['samples']
    samples = " ".join(s)
    return samples

rule pipeline:
    params:
        dimension=config['dimension'],
        samples=_get_samples
    input:
        expand("%s/{barcodes}.fasta" % (DEMUX_DIR),barcodes=BARCODES),
    output:
        "%s/%s_good.fasta" % (BUILD_DIR, PREFIX),
        "%s/%s_partial.fasta" % (BUILD_DIR, PREFIX),
        "%s/%s_poor.fasta" % (BUILD_DIR, PREFIX),
    conda:
        "envs/anaconda.pipeline-env.yaml"
    shell:
        "python pipeline/scripts/pipeline.py --samples {params.samples} --dimension {params.dimension}"
