#!/usr/bin/env python
import os
import argparse

parser = argparse.ArgumentParser(description='Moves all files from numbered subdirectories of workspace/ into workspace, then builds one subdirectory called demux/ in which demultiplexing can take place.')
parser.add_argument('--library', type = str, help="Library to prep")
args = parser.parse_args()

libraries = { '1' : 'usvi-library1-2016-12-10',
              '3' : 'usvi-library3-2017-02-02',
              '4' : 'usvi-library4-2017-03-03',
              '5' : 'usvi-library5-2017-03-14',
              '6' : 'usvi-library6-2017-03-22',
              '7' : 'usvi-library7-1d-2017-03-24',
              '8' : 'usvi-library8-1d-2017-03-31',
              '1_test' : 'usvi-library1-2016-12-10/test/workspace',
              '1-2_test' : 'usvi-library1-2016-12-10/test2/workspace',
              '7_test' : 'usvi-library7-1d-2017-03-24/test/workspace',
              '8_test' : 'usvi-library8-1d-2017-03-31/test/basecalled_reads/workspace' }

assert args.library in libraries.keys(), "%s not a valid library"%(args.library)

# Check for test directories: 1 and 7
if args.library[-4:] != 'test':
    b = '/fh/fast/bedford_t/zika-seq/data/%s/alba121/workspace'%(libraries[args.library])
else:
    b = '/fh/fast/bedford_t/zika-seq/data/%s' % (libraries[args.library])

# Iterate over numbered folders in basecalled directory
for nf in os.listdir(b):
    nfd = '%s/%s'%(b, nf)
    if os.path.isdir(nfd):
        # Move folders out of numbered folder into workspace directory
        for f in os.listdir(nfd):
            of = '%s/%s'%(nfd, f)
            nf = '%s/%s'%(b, f)
            print('mv %s %s'%(of, nf))
            os.rename(of, nf)

# Delete numbered folders
# Note that this will also delete any other subdirectories of workspace
for d in os.listdir(b):
    nfd = '%s/%s'%(b, d)
    if os.path.isdir(nfd):
        os.rmdir(nfd)

# Construct new demux directory
demux = '%s/demux/'%(b)
os.mkdir(demux)
