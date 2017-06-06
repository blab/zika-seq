#!/usr/bin/env python
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--library', type = str, help="Library to prep")
args = parser.parse_args()

libraries = { '7' : 'usvi-library7-1d-2017-03-24' }

assert args.library in libraries.keys(), "%s not a valid library"%(args.library)

b = '/fh/fast/bedford_t/zika-seq/data/%s/test/workspace'%(libraries[args.library])
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
