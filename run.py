#!/usr/bin/env python
import os, subprocess, shutil

libraries = ['usvi-library1-2016-12-10']
barcodes = ['NB03']

for library in libraries:
    for barcode in barcodes:

        # poretools fast5 to fastq
        input_dir = 'data/libraries/' + library + '/basecalled_reads/pass_demultiplex/' + barcode
        if not os.path.exists(input_dir):
            raise NotADirectoryError('Input directory', input_dir, 'does not exist')
        output_dir = 'data/libraries/' + library + '/fastq_reads/'
        output_file = output_dir + barcode + '.fastq'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        f = open(output_file, "w")
        call = map(str, ['poretools', 'fastq', '--type', '2D', input_dir])
        print(" ".join(call))
        subprocess.call(call, stdout=f)

        # marginAlign to reference
        input_file = 'data/libraries/' + library + '/fastq_reads/' + barcode + '.fastq'
        output_dir = 'data/libraries/' + library + '/sam_reads/'
        output_file = output_dir + barcode + '.sam'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        # marginAlign input.fastq ref.fasta out.sam
        call = map(str, ['marginAlign', input_file, 'refs/Zika_FP.fasta', output_file])
        if os.path.exists('jobTree'):
            shutil.rmtree('jobTree')
        print(" ".join(call))
        subprocess.call(call)
        if os.path.exists('jobTree'):
            shutil.rmtree('jobTree')

        # samtools convert sam to sorted bam
        input_file = 'data/libraries/' + library + '/sam_reads/' + barcode + '.sam'
        output_dir = 'data/libraries/' + library + '/bam_reads/'
        output_file = output_dir + barcode + '.bam'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        # samtools view -bS file.sam | samtools sort -o file_sorted
        call = 'samtools view -bS ' + input_file + ' | samtools sort -o ' + output_file
        print(call)
        subprocess.call(call, shell=True)

        # samtools create index from sorted bam
        input_file = 'data/libraries/' + library + '/bam_reads/' + barcode + '.bam'
        output_file = 'data/libraries/' + library + '/bam_reads/' + barcode + '.bai'
        # samtools index test_sorted.bam test_sorted.bai
        call = map(str, ['samtools', 'index', input_file, output_file])
        print(" ".join(call))
        subprocess.call(call)
