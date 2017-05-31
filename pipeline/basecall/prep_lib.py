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
for nf in os.listdir(b):
    nfd = '%s/%s'%(b, nf)
    if os.path.isdir(nfd):
        for f in os.listdir(nfd):
            of = '%s/%s'%(nfd, f)
            nf = '%s/%s'%(b, f)
            print('mv %s %s'%(of, nf))
            os.rename(of, nf)
for d in os.listdir(b):
    nfd = '%s/%s'%(b, d)
    if os.path.isdir(nfd):
        os.rmdir(nfd)
demux = '%s/demux/'%(b)
os.mkdir(demux)
