#!/usr/bin/env python
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--library', type = str, help="Library to prep")
args = parser.parse_args()

libraries = { '3' : 'usvi-library3-2017-02-02',
              '4' : 'usvi-library4-2017-03-03',
              '5' : 'usvi-library5-2017-03-14',
              '6' : 'usvi-library6-2017-03-22',
              '8' : 'usvi-library8-1d-2017-03-31' }

assert args.library in libraries.keys(), "%s not a valid library"%(args.library)

b = '/fh/fast/bedford_t/zika-seq/data/%s/basecalled_reads/workspace'%(libraries[args.library])
for bc in os.listdir(b):
    bcd = '%s/%s'%(b, bc)
    if bcd[-6:] != '.fast5':
        for f in os.listdir('%s/0/'%(bcd)):
            if f[-6:] == '.fast5':
                of = '%s/0/%s'%(bcd, f)
                nf = '%s/%s'%(bcd, f)
                print('mv %s %s'%(of, nf))
                os.rename(of, nf)
