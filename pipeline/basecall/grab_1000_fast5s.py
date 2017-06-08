#!/usr/bin/env python
import os, shutil

from_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/basecalled_reads/'
to_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/basecalled_reads/workspace/demux/'
count = 0

for f in os.listdir(from_dir):
    if f[-6:] == '.fast5' and count < 1000:
        old = from_dir + f
        new = to_dir + f
        print()
        print('%s: mv %s %s'%(count, old, new))
        # shutil.copyfile()
        count += 1
