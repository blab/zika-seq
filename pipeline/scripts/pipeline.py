#!/usr/bin/env python
import argparse, csv, subprocess, time
import sys
from Bio import SeqIO
import os

def sample_to_run_data_mapping(samples_dir):
    '''
    return dict
    each key is string "sample_id"
    each value is a list of tuples ("library", "barcode")
    '''
    runs_file = samples_dir + "runs.tsv"
    sr_mapping = {}
    with open(runs_file) as tsv:
        for row in csv.DictReader(tsv, delimiter="\t"):
            sample = row["sample_id"]
            rb_pair = (row["run_name"], row["barcode_id"])
            if sample not in sr_mapping:
                sr_mapping[sample] = []
            sr_mapping[sample].append(rb_pair)
    return sr_mapping

def sample_to_metadata_mapping(samples_dir):
    '''
    return dict
    each key is string "sample_id"
    each value is a list of metadata ordered as
    ["strain", "sample_id", "collect_date", "country", "division", "location"]
    '''
    metadata_file = samples_dir + "samples.tsv"
    sm_mapping = {}
    with open(metadata_file) as tsv:
        for row in csv.DictReader(tsv, delimiter="\t"):
            sample = row["sample_id"]
            metadata = [row["strain"], row["sample_id"], row["collection_date"],
                row["country"], row["division"], row["location"]]
            sm_mapping[sample] = metadata
    return sm_mapping

def construct_sample_fastas(sr_mapping, data_dir, build_dir):
    ''' Construct one fasta for each sample from two input fastas.
    '''
    import gzip
    for sample in sr_mapping:
        # Grab a matched pair of barcode fastas; global paths
        fastas = [ '%s%s/alba121/workspace/demux/%s.fasta' % (data_dir, run, barcode) for (run, barcode) in sr_mapping[sample] ]
        for fasta in fastas:
            if fasta.endswith('na.fasta'):
                fastas.remove(fasta)
                print('Remvoed 1 fasta ending in na.fasta')
        print(fastas)
        assert len(fastas) == 2, 'Expected 2 .fasta files for %s, instead found %s.\nCheck that they are present and gzipped in %s%s/basecalled_reads/workspace/demux/' % (sample, len(fastas), data_dir, sr_mapping[sample][0])
        complete_fasta = '%s%s_complete.fasta' % (build_dir, sample)
        with open(complete_fasta, 'w+') as f:
            with open(fastas[0], 'r') as f1:
                print('Writing %s to %s' % (fastas[0],complete_fasta))
                content = f1.read()
                f.write(content)
            with open(fastas[1], 'r') as f2:
                print('Writing %s to %s' % (fastas[1],complete_fasta))
                content = f2.read()
                f.write(content)

        final_fasta = '%s%s.fasta' % (build_dir, sample)
        with open(final_fasta, 'w+') as f:
            (run, barcode) = sr_mapping[sample][0]
            sed_str = '%s%s/alba121/workspace' % (data_dir, run)
            sed_str = sed_str.split('/')
            sed_str = '\/'.join(sed_str)
            call = 'sed \'s/\.\./%s/\' %s%s_complete.fasta' % (sed_str, build_dir, sample)
            print(" > ".join([call, final_fasta]))
            subprocess.call(call, shell=True, stdout=f)

def process_sample_fastas(sm_mapping, build_dir, dimension):
    ''' Run fasta_to_consensus script to construct consensus files.
    TODO: Make sure that this runs after changes to inputs and fasta_to_consensus on 1d reads
    '''
    for sample in sm_mapping:
        print("* Processing " + sample)
        # build consensus
        sample_stem = build_dir + sample
        if dimension == '2d':
            call = ['pipeline/scripts/fasta_to_consensus_2d.sh', '/fh/fast/bedford_t/zika-seq/pipeline/refs/KJ776791.2.fasta', sample_stem, '/fh/fast/bedford_t/zika-seq/pipeline/metadata/v2_500.amplicons.ver2.bed']
        elif dimension == '1d':
            call = ['pipeline/scripts/fasta_to_consensus_1d.sh', '/fh/fast/bedford_t/zika-seq/pipeline/refs/KJ776791.2.fasta', sample_stem, '/fh/fast/bedford_t/zika-seq/pipeline/metadata/v2_500.amplicons.ver2.bed']
        print(" ".join(call))
        subprocess.call(call)
        # annotate consensus
        # >ZBRD116|ZBRD116|2015-08-28|brazil|alagoas|arapiraca|minion
        print('#############\n')
        fasta_header = ">" + "|".join(sm_mapping[sample])
        fasta_header += "|minion"
        replacement = r"\~^>~s~.*~" + fasta_header + "~" # ~ rather than / to avoid conflict with strain names
        input_file = build_dir + sample + ".consensus.fasta"
        output_file = "temp.fasta"
        f = open(output_file, "w")
        call = ['sed', replacement, input_file]
        print(" ".join(call) + " > " + output_file)
        subprocess.call(call, stdout=f)
        call = ['mv', output_file, input_file]
        print(" ".join(call))
        subprocess.call(call)
        print("")

def gather_consensus_fastas(sm_mapping, build_dir, prefix):
    '''
    Gather consensus files into genomes with 'partial' (50-80 percent coverage)
    and good (>80 percent coverage) coverage
    '''
    # identify partial and good samples
    print("* Concatenating consensus fastas")
    partial_samples = []
    good_samples = []
    poor_samples = []
    for sample in sm_mapping:
        consensus_file = build_dir + sample + ".consensus.fasta"
        with open(consensus_file) as f:
            lines = f.readlines()
        seq = lines[1]
        coverage = 1 - seq.count("N") / float(len(seq))
        print(seq.count("N")) #DEBUG
        print(len(seq)) #DEBUG
        print("COVERAGE: "+ str(coverage)) #DEBUG
        if coverage >= 0.5 and coverage < 0.8:
            partial_samples.append(sample)
        elif coverage >= 0.8:
            good_samples.append(sample)
        else:
            poor_samples.append(sample)
    # sort samples
    partial_samples.sort()
    good_samples.sort()
    poor_samples.sort()
    print("Good samples: " + " ".join(good_samples))
    print("Partial samples: " + " ".join(partial_samples))
    print("Poor samples: " + " ".join(poor_samples))
    # concatenate partial samples
    input_file_list = [build_dir + sample + ".consensus.fasta" for sample in partial_samples]
    output_file = build_dir + prefix + "_partial.fasta"
    f = open(output_file, "w")
    call = ['cat'] + input_file_list
    print(" ".join(call) + " > " + output_file)
    if len(input_file_list) >= 1:
        subprocess.call(call, stdout=f)
    # concatenate good samples
    input_file_list = [build_dir + sample + ".consensus.fasta" for sample in good_samples]
    output_file = build_dir + prefix + "_good.fasta"
    f = open(output_file, "w")
    call = ['cat'] + input_file_list
    print(" ".join(call) + " > " + output_file)
    subprocess.call(call, stdout=f)
    # concatenate poor samples
    print("Poor samples: " + " ".join(good_samples))
    input_file_list = [build_dir + sample + ".consensus.fasta" for sample in poor_samples]
    output_file = build_dir + prefix + "_poor.fasta"
    f = open(output_file, "w")
    call = ['cat'] + input_file_list
    print(" ".join(call) + " > " + output_file)
    subprocess.call(call, stdout=f)
    print("")

def overlap(sr_mapping, build_dir):
    # prepare sorted bam files for coverage plots
    for sample in sr_mapping:
        # samtools depth <name.sorted.bam> > <name.coverage>
        bamfile = build_dir + sample + '.sorted.bam'
        coveragefile = build_dir + sample + '.coverage'
        with open(coveragefile, 'w+') as f:
            call = ['samtools', 'depth', bamfile]
            print(" ".join(call + ['>', coveragefile]))
            subprocess.call(call, stdout=f)

        chfile = build_dir + sample + '.chr1.coverage'
        with open(coveragefile, 'r') as f:
            f.readline()
            line = f.readline()
            l = line.split('\t')
            chromosome_name = l[0]
        call = "awk '$1 == \"" + chromosome_name + "\" {print $0}' " + coveragefile + " > " + chfile
        print(call)
        subprocess.call([call], shell=True)
        call = "Rscript /fh/fast/bedford_t/zika-seq/pipeline/scripts/depth_coverage.R --inFile " + chfile + " --outPath " + build_dir + " --name " + sample
        print(call)
        subprocess.call([call], shell=True)
        print("")

def per_base_error_rate(sr_mapping, build_dir):
    '''
    Calculate per-base error rates by walking through VCF files for each sample.
    TODO: make this work
    '''
    length = 10794.0 # TODO: Make sure this always works or is variable
    for sample in sr_mapping:
        error = 0
        vcf = build_dir + sample + '.vcf'
        if vcf in os.listdir(build_dir):
            with open(vcf) as f:
                lines = f.readlines()
                if len(lines) > 1:
                    for line in lines:
                        l = line.split('\t')
                        alt = len(l[4])
                        error += alt
            outfile = build_dir + sample + '.error'
            error = error / length
            with open(outfile, 'w+') as f:
                f.write('Error rate: ' + str(error))

if __name__=="__main__":
    parser = argparse.ArgumentParser( description = "Bioinformatic pipeline for generating consensus genomes from demultiplexed Zika fastas" )
    parser.add_argument( '--data_dir', type = str, default = "/fh/fast/bedford_t/zika-seq/data/",
                            help="directory containing data; default is \'/fh/fast/bedford_t/zika-seq/data/\'")
    parser.add_argument( '--samples_dir', type = str, default = "/fh/fast/bedford_t/zika-seq/samples/",
                            help="directory containing samples and runs TSVs; default is \'/fh/fast/bedford_t/zika-seq/samples/\'" )
    parser.add_argument( '--build_dir', type = str, default = "/fh/fast/bedford_t/zika-seq/build/",
                            help="directory for output data; default is \'/fh/fast/bedford_t/zika-seq/build/\'" )
    parser.add_argument('--prefix', type = str, default = "ZIKA_USVI",
                            help="string to be prepended onto all output consensus genome files; default is \'ZIKA_USVI\'")
    parser.add_argument('--samples', type = str, default = None, nargs='*',
                            help="sample to be run")
    parser.add_argument('--dimension', type = str, default = '2d',
                            help="dimension of library to be fun; options are \'1d\' or \'2d\', default is \'2d\'")
    parser.add_argument('--run_steps', type = int, default = None, nargs='*',
                            help="Numbered steps that should be run (i.e. 1 2 3):\n\t1. Construct sample fastas \n\t2. Process sample fastas \n\t3. Gather consensus fastas \n\t 4. Generate overlap graphs \n\t5. Calculate per-base error rates")
    params = parser.parse_args()

    assert params.dimension in [ '1d', '2d' ], "Unknown library dimension: options are \'1d\' or \'2d\'."

    dd = params.data_dir
    bd = params.build_dir
    sd = params.samples_dir

    # This will reduce my own frustration
    if dd[-1] != '/':
        dd += '/'
    if bd[-1] != '/':
        bd += '/'
    if sd[-1] != '/':
        sd += '/'

    start_time = time.time()

    # Parse samples directory and remove unwanted samples
    # TODO: This can be run more efficiently in the future by making sample_to_run_data_mapping and sample_to_metadata_mapping smarter
    sr_mapping = sample_to_run_data_mapping(sd)
    sm_mapping = sample_to_metadata_mapping(sd)
    tmp_sr = { s: sr_mapping[s] for s in params.samples }
    tmp_sm = { s: sm_mapping[s] for s in params.samples }
    sr_mapping = tmp_sr
    sm_mapping = tmp_sm

    # Only run specified run_steps.
    # If pipeline is modified, add helper function here and put its name in pipeline below
    def csf():
        construct_sample_fastas(sr_mapping, dd, bd)
    def psf():
        process_sample_fastas(sm_mapping, bd, params.dimension)
    def gcf():
        gather_consensus_fastas(sm_mapping, bd, params.prefix)
    def go():
        overlap(sm_mappng, bd)
    def pber():
        per_base_error_rate(sr_mapping, bd)

    if params.run_steps is not None:
        for index in params.run_steps:
            assert index in [1,2,3,4,5], 'Unknown step number %s, options are 1, 2, 3, 4, or 5.' % (index)
    pipeline = [ csf, psf, gcf, go, pber ]
    if params.run_steps is None:
        for fxn in pipeline:
            fxn()
    else:
        for index in params.run_steps:
            pipeline[index-1]()

    time_elapsed = time.time() - start_time
    m, s = divmod(time_elapsed, 60)
    h, m = divmod(m, 60)
    print('Total runtime: %d:%02d:%02d' % (h, m, s))
