#!/usr/bin/env python
import shutil, argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Completely remove the contents of a library\'s basecalled workspace directory.')
    parser.add_argument('--library', default=None, type=str, help="library that needs a clean workspace")

    args = parser.parse_args()

    datadirs = {    '1' : '/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/basecalled_reads/workspace/',
                    '3' : '/fh/fast/bedford_t/zika-seq/data/usvi-library3-2017-02-02/basecalled_reads/workspace/',
                    '4' : '/fh/fast/bedford_t/zika-seq/data/usvi-library4-2017-03-03/basecalled_reads/workspace/',
                    '5' : '/fh/fast/bedford_t/zika-seq/data/usvi-library5-2017-03-14/basecalled_reads/workspace/',
                    '6' : '/fh/fast/bedford_t/zika-seq/data/usvi-library6-2017-03-22/basecalled_reads/workspace/',
                    '7' : '/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/basecalled_reads/workspace/',
                    '8' : '/fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/basecalled_reads/workspace/'
                }

    yes = [ 'yes', 'y' ]

    assert args.library in datadirs.keys(), "Unknown library."

    remove_directory = datadirs[args.library]

    check = raw_input('Are you sure you want to delete %s? (y/n) '%(remove_directory))

    if check.lower() in yes:
        print('Removing %s'%(remove_directory))
        shutil.rmtree(remove_directory)
    else:
        print('Not deleting any directory.')
