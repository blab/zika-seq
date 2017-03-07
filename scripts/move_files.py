import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', help='directory that needs to be sorted')
parser.add_argument('-o','--output', default=None, help='directory to put raw reads')

def move_files(indir, outdir):
    rr = indir
    for number in os.listdir(rr):
        print('Emptying folder ' + str(number))
        nrr = rr + str(number) + '/'
        print(nrr)
        for f in os.listdir(nrr):
            oldf = nrr + f
            if outdir:
                newf = outdir + f
            else:
                newf = rr + f
            print('Moving ' + oldf + ' to ' + newf)
            os.rename(oldf, newf)
        os.rmdir(nrr)

if __name__=="__main__":
    args = parser.parse_args()
    move_files(args.input, args.output)
