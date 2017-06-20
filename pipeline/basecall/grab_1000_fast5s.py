#!/usr/bin/env python
import os, shutil

parser = argparse.ArgumentParser(description='Move 1000 fast5 files from either: \n\tlibrary1/basecalled_reads/workspace/ to library1/basecalled_reads/workspace/demux/poretools_test\n\tlibrary7/basecalled_reads/workspace/ to library7/basecalled_reads/workspace/demux/poretools_test/')
parser.add_argument('--dimension', type = str, help="Dimension of reads to grab: 1d or 2d")
args = parser.parse_args()

if args.dimension == '2d':
    from_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/basecalled_reads/workspace/'
    to_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/basecalled_reads/workspace/demux/poretools_test/'
else:
    from_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/test/workspace/'
    to_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/test/workspace/demux/poretools_test/'

count = 0
for f in os.listdir(from_dir):
    # Only move 1000 files
    while count < 1000:
        if f[-6:] == '.fast5':
            old = from_dir + f
            new = to_dir + f
            print('%s: mv %s %s'%(count, old, new))
            shutil.copyfile(old, new)
            count += 1
