#!/usr/bin/env python
import os

b = '/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/basecalled_reads/workspace'
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
demux = '%s/demux/'
os.mkdir(demux)
