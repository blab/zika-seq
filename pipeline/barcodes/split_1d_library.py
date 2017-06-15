#!/bin/bash/python
'''
split_1d_library.py
'''
import argparse, os

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Script to move fast5 files into batches of 500,000, within numbered folders.')
    parser.add_argument('--libraryPath', type=str, help='Path to the library that needs to be split up.')
    parser.add_argument('--run', default=False, action='store_true', help='Run the move command')
    parser.add_argument('--verbosity', default=1, type=int, help='Level of verbosity in output: 0, 1, or 2')
    args = parser.parse_args()
    # Default number if files/directory is 500,000
    dirsize = 500000
    count = 0
    lib = args.libraryPath

    assert args.verbosity in [0, 1, 2], "Unknown verbosity level. Options are 0, 1, or 2."

    # Only look at fast5 files in the library
    for f in os.listdir(lib):
        if f[-6:] == '.fast5':
            # Make new directories once 500,000 reads have gone into the old directory
            dirnum = str(count // dirsize)
            dirname = lib + dirnum + '/'
            if count % dirsize == 0:
                if not os.path.exists(dirname):
                    if args.verbosity in [1, 2]:
                        print('mkdir %s' % (dirname))
                    os.makedirs(dirname)
                else:
                    if args.verbosity in [1, 2]:
                        print('%s already exists.' % (dirname))
            # Once a new directory is made, move files
            old = lib + f
            new = dirname + f
            if args.verbosity == 2:
                print('%s: mv %s %s' % (count, old, new))
            # Argument for debug
            if args.run:
                os.rename(old, new)

            count += 1
