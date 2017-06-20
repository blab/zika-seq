#!/usr/bin/env python
import os
'''
Count how many reads were demultiplexeded for each library using Albacore.
Note that Albacore native demuxing is no longer being used.
'''


out_file = '/fh/fast/bedford_t/zika-seq/pipeline/basecall/ratios.txt'
prefix = '/fh/fast/bedford_t/zika-seq/data/'

# Names of Albacore output directories
barcodes = [ 'barcode01', 'barcode02', 'barcode03', 'barcode04',
             'barcode05', 'barcode06', 'barcode07', 'barcode08',
             'barcode09', 'barcode10', 'barcode11', 'barcode12' ]

# Our library data paths
libraries = { '1' : 'usvi-library1-2016-12-10/basecalled_reads/workspace',
              '3' : 'usvi-library3-2017-02-02/basecalled_reads/workspace',
              '4' : 'usvi-library4-2017-03-03/basecalled_reads/workspace',
              '5' : 'usvi-library5-2017-03-14/basecalled_reads/workspace',
              '6' : 'usvi-library6-2017-03-22/basecalled_reads/workspace',
              '7' : 'usvi-library7-1d-2017-03-24/basecalled_reads/workspace',
              '8' : 'usvi-library8-1d-2017-03-31/basecalled_reads/workspace' }

with open(out_file, 'w+') as f:
    f.write('Library\tBasecalled\tUnclassified\tRatio\n')
    for library in libraries.keys():
        barcode_count = float(sum( [ len(os.listdir('%s%s/%s/'%(prefix, libraries[library], barcode))) for barcode in barcodes ] ))
        unclassified_count = float(len(os.listdir('%s%s/unclassified'%(prefix, libraries[library]))))
        ratio = barcode_count / unclassified_count
        f.write('%s\t%s\t%s\t%s'%(library, str(barcode_count), str(unclassified_count), str(ratio)))
