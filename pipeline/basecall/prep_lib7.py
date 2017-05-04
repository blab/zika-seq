#!/usr/bin/env python
import os

b = '/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/basecalled_reads/workspace/'
for bc in os.listdir(b):
    bcd = '%s/%s'%(b, bc)
    for f in os.listdir('%s/0/'%(bcd)):
        of = '%s/0/%s'%(bcd, f)
        nf = '%s/%s'%(bcd, f)
        print('mv %s %s'%(of, nf))
        # os.rename(of, nf)
