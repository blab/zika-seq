#!/usr/bin/env python
import os, shutil, sys

from_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/test/workspace/'
to_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/test/workspace/demux/poretools_test/'
count = 0

for f in os.listdir(from_dir):
    if f[-6:] == '.fast5' and count < 1000:
        old = from_dir + f
        new = to_dir + f
        print('%s: mv %s %s'%(count, old, new))
        shutil.copyfile(old, new)
        count += 1
    else:
        sys.exit()
