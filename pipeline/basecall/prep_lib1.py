#!/usr/bin/env python
import os

b = '/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/basecalled_reads/workspace'
for bc in os.listdir(b):
    bcd = '%s/%s'%(b, bc)
    if os.path.isdir(bcd):
        for f in os.listdir('%s/0/'%(bcd)):
            of = '%s/0/%s'%(bcd, f)
            nf = '%s/%s'%(bcd, f)
            print('mv %s %s'%(of, nf))
            os.rename(of, nf)
