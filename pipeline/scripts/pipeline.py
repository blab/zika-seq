#!/usr/bin/env python
import argparse, csv, subprocess, time
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

def construct_sample_fastas(sr_mapping, data_dir, build_dir, logfile, dimension):
    '''
    Use nanopolish to construct a single fasta for all reads from a sample
    '''
    # runs = set()
    # for sample in sr_mapping:
    #     for (run, barcode) in sr_mapping[sample]:
    #         runs.add(run)
    # for r in runs:
    #     fail_folder = data_dir + r + "/basecalled_reads/fail/"
    #     dmname = r + '_fail_demultiplex.fasta'
    #     demultiplex_file = build_dir + dmname
    #     print("Demultiplexing: "+demultiplex_file)
    #     if dmname not in os.listdir(build_dir):
    #         f = open(demultiplex_file, "w+")
    #         # Take this out if we ignore fail reads
    #         call = [ 'poretools', 'fasta', '--type', '2D', fail_folder ]
    #         print(" ".join(call))
    #         subprocess.call(call,stdout=f)
    #         f.close()
    #         call = [ 'barcodes/demultiplex.py', '--barcodes', '/zibra/zika-pipeline/barcodes/barcodes.fasta', demultiplex_file]
    #         print(" ".join(call))
    #         subprocess.call(call)
    #         with open(logfile, 'a') as f:
    #             f.write(time.strftime('[%H:%M:%S] Done demultiplexing ' + run + '\n'))
    #
    # # clean up the demultiplexed files.
    # for fl in os.listdir(build_dir):
    #     if fl[0] == '>':
    #         print("Fixing name of " + fl)
    #         fname = fl[1:]
    #         os.rename(build_dir + fl, build_dir + fname)
    #         with open(logfile, 'a') as f:
    #             f.write(time.strftime('[%H:%M:%S] Fixed name of ' + fname + '\n'))

    import glob
    for sample in sr_mapping:
        bc = []
        fail_file = build_dir + sample + '_fail_demultiplex.fasta'
        f = open(fail_file, "w+")
        for (run, barcode) in sr_mapping[sample]:
            bc += glob.glob(build_dir + barcode + '_' + run + '_fail_demultiplex.fasta')
        call = ['cat'] + bc
        if len(bc) >= 1:
            print(" ".join(call + ['>', fail_file]))
            subprocess.call(call, stdout=f)
        else:
            with open(logfile, 'a') as f:
                f.write("Unable to cat, no demultiplexed files for " + sample)
        with open(logfile, 'a') as f:
            f.write(time.strftime('[%H:%M:%S] Done writing fail fasta for ' + sample + '\n'))

    # Pass reads
    for sample in sr_mapping:
        print("* Extracting " + sample)
        # nanopolish extract each run/barcode pair
        for (run, barcode) in sr_mapping[sample]:
            input_dir = data_dir + run + "/pass/" + barcode # Update this to /basecalled_reads/workspace/
            output_file = build_dir + sample + "_" + run + "_" + barcode + ".fasta"
            f = open(output_file, "w")
            if output_file not in os.listdir(build_dir):
                if dimension == '2d':
                    call = ['nanopolish', 'extract', '--type', '2d', input_dir]
                elif dimension == '1d':
                    call = ['nanopolish', 'extract', '--type', '2d', input_dir]
                print(" ".join(call) + " > " + output_file)
                subprocess.call(call, stdout=f)
            else:
                with open(logfile, 'a') as f:
                    f.write(time.strftime('[%H:%M:%S] ' + output_file  + ' already in ' + build_dir + '\n'))

        # concatenate to single sample fasta
        input_file_list = [build_dir + sample + "_" + run + "_" + barcode + ".fasta"
            for (run, barcode) in sr_mapping[sample]]
        # add demultiplex file
        failed = build_dir + sample + '_fail_demultiplex.fasta'
        input_file_list.append(failed)
        output_file = build_dir + sample + ".fasta"
        f = open(output_file, "w")
        call = ['cat'] + input_file_list
        # Check if any files exist to prevent freezing
        if len(input_file_list) >= 1:
            print(" ".join(call) + " > " + output_file )
            subprocess.call(call, stdout=f)
            with open(logfile, 'a') as f:
                f.write(time.strftime('[%H:%M:%S] Done writing complete fasta for ' + sample + '\n'))
        else:
            with open(logfile, 'a') as f:
                f.write("Unable to cat, no fasta files available for " + sample)
        print("")

def process_sample_fastas(sm_mapping, build_dir, logfile, dimension):
    '''
    Run fasta_to_consensus script to construct consensus files
    '''
    for sample in sm_mapping:
        print("* Processing " + sample)
        # build consensus
        sample_stem = build_dir + sample
        if dimension == '2d':
            call = ['pipeline/scripts/fasta_to_consensus_2d.sh', 'pipeline/refs/KJ776791.2.fasta', sample_stem, 'pipeline/metadata/v2_500.amplicons.ver2.bed']
        elif dimension == '1d':
            call = ['pipeline/scripts/fasta_to_consensus_1d.sh', 'pipeline/refs/KJ776791.2.fasta', sample_stem, 'pipeline/metadata/v2_500.amplicons.ver2.bed']
        print(" ".join(call))
        subprocess.call(call)
        # annotate consensus
        # >ZBRD116|ZBRD116|2015-08-28|brazil|alagoas|arapiraca|minion
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
        with open(logfile, 'a') as f:
            f.write(time.strftime('[%H:%M:%S] Consensus fasta completed for ' + sample + '\n'))
        print("")

def gather_consensus_fastas(sm_mapping, build_dir, prefix, logfile):
    '''
    Gather consensus files into genomes with 'partial' (50-80% coverage)
    and good (>80% coverage) coverage
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
    with open(logfile, 'a') as f:
        f.write(time.strftime('[%H:%M:%S] Done gathering consensus fastas\n'))
    print("")

def overlap(sr_mapping, build_dir, logfile):
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
        call = "Rscript scripts/depth_coverage.R -i " + chfile + " -o " + build_dir + " -n " + sample
        print(call)
        subprocess.call([call], shell=True)
        print("")
        with open(logfile, 'a') as f:
            f.write(time.strftime('[%H:%M:%S] Done drawing overlap graphs for ' + sample + '\n'))

def per_base_error_rate(sr_mapping, build_dir, logfile):
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
    with open(logfile, 'a') as f:
        f.write(time.strftime('[%H:%M:%S] Done calculating per-base error rates ' + sample + '\n'))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "process data")
    parser.add_argument('--data_dir', type = str, default = "data/")
    parser.add_argument('--samples_dir', type = str, default = "samples/")
    parser.add_argument('--build_dir', type = str, default = "build/")
    parser.add_argument('--prefix', type = str, default = "ZIKA_USVI")
    parser.add_argument('--samples', type = str, nargs='*', default = None)
    parser.add_argument('--dimension', type = str, default = '2d')
    params = parser.parse_args()

    assert params.dimension in [ '1d', '2d' ], "Unknown library dimension"

    logfile = params.build_dir + 'log.txt'
    start_time = time.time()
    with open(logfile,'w+') as f:
        f.write(time.strftime('Pipeline started on %Y-%m-%d at %H:%M:%S\n'))

    sr_mapping = sample_to_run_data_mapping(params.samples_dir)
    sm_mapping = sample_to_metadata_mapping(params.samples_dir)
    if params.samples:
        for sample in sr_mapping.keys():
            if sample not in params.samples:
                sr_mapping.pop(sample, None)
            if sample not in params.samples:
                sm_mapping.pop(sample, None)

    with open(logfile,'a') as f:
        f.write('Samples:\n')
        for sample in sr_mapping:
            f.write(sample+'\n')
            for (run, barcode) in sr_mapping[sample]:
                f.write('\t('+run+', '+barcode+')\n')

    construct_sample_fastas(sr_mapping, params.data_dir, params.build_dir, logfile, params.dimension)
    process_sample_fastas(sm_mapping, params.build_dir, logfile, params.dimension)
    gather_consensus_fastas(sm_mapping, params.build_dir, params.prefix, logfile)
    overlap(sm_mapping, params.build_dir, logfile)
    per_base_error_rate(sr_mapping, params.build_dir, logfile)

    time_elapsed = time.time() - start_time
    m, s = divmod(time_elapsed, 60)
    h, m = divmod(m, 60)
    with open(logfile,'a') as f:
        f.write(time.strftime('Pipeline completed on %Y-%m-%d at %H:%M:%S\n'))
        f.write('Total runtime: %d:%02d:%02d' % (h, m, s))
