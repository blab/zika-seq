#!/usr/bin/env python
import os

lib3dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library3-2017-02-02/'
lib3raw = lib3dir + 'raw_reads'

for f in os.listdir(lib3dir):
    if f[-6:] == '.fast5':
        old = lib3dir + f
        new = lib3raw + f
        print('mv %s %s'%(old, new))
        os.rename(old, new)
