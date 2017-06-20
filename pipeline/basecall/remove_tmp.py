#!/usr/bin/env python
import argparse, os

parser = argparse.ArgumentParser(description='Remove .tmp prefix of library 8 files that were not completely transferred at minKNOW runtime.')
parser.add_argument('-i', '--input', default='/fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/tmp/tmp/')
parser.add_argument('-o', '--output', default='/fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/raw_reads/')

def rename_contents(inpath, outpath):
    for old_fname in os.listdir(inpath):
        if old_fname[-4:] == '.tmp':
            new_fname = old_fname[:-4]
            old_path = inpath + old_fname
            new_path = outpath + new_fname
            try:
                os.rename(old_path, new_path)
            except:
                pass

def make_directory(path, folder):
    directory = path + folder
    if not os.path.exists(directory):
        os.makedirs(directory)

if __name__=="__main__":
    args = parser.parse_args()

    for folder in os.listdir(args.input):
        make_directory(args.output, folder)
        infolder = args.input + folder + '/'
        outfolder = args.output + folder + '/'
        rename_contents(infolder, outfolder)
