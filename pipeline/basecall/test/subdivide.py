#!/usr/bin/env python
import os, argparse


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--library', type=str, help='input library that needs to be subdivided')
    parser.add_arguemnt('--real', default=False, action='store_true', help='actually run the command, not just print the commands')

    args = parser.parse_args()

    lib = args.library
    print(lib)

    dircount = 0
    count = 0

    starter = lib + '0/'
    if args.real:
        if not os.path.exists(starter):
            os.makedirs(starter)

    for f in os.listdir(lib):
        # Make sure that we are only iterating over .fast5 files
        if f[-6:] == '.fast5':
            if count <= 3999:
                oldpath = '%s%s'%(lib, f)
                newpath = '%s%s/%s'%(lib, str(dircount), f)
                if (count % 1000) == 0:
                    print('%s: mv %s %s'%(count, oldpath, newpath))
                if args.real:
                    os.rename(oldpath, newpath)
                count += 1
            else:
                dircount += 1
                newdir = '%s%s/'%(lib, str(dircount))
                print('Count = %s; making new directory: %s.'%(count, newdir))
                if args.real:
                    if not os.path.exists(directory):
                        os.makedirs(directory)
                oldpath = '%s%s'%(lib, f)
                newpath = '%s%s/%s'%(lib, str(dircount), f)
                print('%s: mv %s %s'%(count, oldpath, newpath))
                if args.real:
                    os.rename(oldpath, newpath)
                count = 0
