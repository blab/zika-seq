#!/usr/bin/env python
import subprocess, os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inpath', default=None--, help="path to input directory containing non-basecalled reads")
parser.add_argument('-o', '--outpath', default=None, help="path to output directory for basecalled reads")
parser.add_argument('--dimension', default='1d', help="dimension of sequenced library; 1d or 2d")
parser.add_argument('--email', default=None, help="email address for sbatch notifications")

parser.add_argument('--test', default=False, action='store_true', help="make sure that albacore is working correctly; do not run from SBATCH")
parser.add_argument('--dryrun', default=False, action='store_true', help="print commands to be run, but do not call them")
args = parser.parse_args()

if args.dimension == '1d':
    w = 'r94_450bps_linear.cfg'
else:
    w = 'r94_250bps_2d.cfg'

if args.test:
    call = 'sbatch --wrap=\"read_fast5_basecaller.py -h\"'
    print(call)
    try:
        subprocess.call(call, shell=True)
    except:
        print("Unable to run Albacore. Ensure that proper requisite module is loaded with:\n\nmodule load Python/3.5.2-foss-2016b-fh1\n")
else:
    assert args.inpath, "No input path given"
    assert args.outpath, "No output path given"
    assert args.email, "No email address given; necessary for SBATCH"

    for dirname in os.listdir(args.inpath):
        acall = 'read_fast5_basecaller.py -i %s/%s -t 6 -c %s -r --barcoding -s %s'%(args.inpath, dirname, w, args.outpath)
        call = [ 'sbatch', '--time=24:00:00', '--mem=10000', '--mail-type=FAIL', '--mail-user=%s'%(args.email), '--wrap="%s"'%(acall) ]
        print(" ".join(call))
        if not args.dryrun:
            subprocess.call(" ".join(call), shell=True)
