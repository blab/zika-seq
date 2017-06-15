#!/bin/bash/python
'''
split_1d_library.py
'''
import argparse, os

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Script to move fast5 files into batches of 500,000, within numbered folders.')
    parser.add_argument('--libraryPath', type=str, help='Path to the library that needs to be split up.')
    parser.add_argument('--run', default=False, action='store_true', help='Run the move command')
    args = parser.parse_args
    dirsize = 500000
    count = 0
    lib = args.libraryPath

    for f in os.listdir(lib):
        if f[-6:] == '.fast5':
            if count % dirsize == 0:
                os.makedir(lib + str(count // dirsize) + '/')
            old = lib + f
            new = lib + str(count // dirsize) + '/' + f
            print('%s: mv %s %s' % (count, old, new))
            if args.run:
                os.rename(old, new)
            count += 1
