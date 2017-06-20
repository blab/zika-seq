#!/usr/bin/env python
import os, argparse


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--library', type=str, help='input library that needs to be subdivided')
    parser.add_argument('--real', default=False, action='store_true', help='actually run the command, not just print the commands')
    args = parser.parse_args()

    libraries = { '1' : 'usvi-library1-2016-12-10',
                  '3' : 'usvi-library3-2017-02-02',
                  '4' : 'usvi-library4-2017-03-03',
                  '5' : 'usvi-library5-2017-03-14',
                  '6' : 'usvi-library6-2017-03-22',
                  '7' : 'usvi-library7-1d-2017-03-24',
                  '8' : 'usvi-library8-1d-2017-03-31',
                  '1_test' : 'usvi-library1-2016-12-10/test/workspace',
                  '7_test' : 'usvi-library7-1d-2017-03-24/test/workspace' }

    assert args.library in libraries.keys(), "%s not a valid library, options are %s"%(args.library, str(libraries.keys()))
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
                    if not os.path.exists(newdir):
                        os.makedirs(newdir)
                oldpath = '%s%s'%(lib, f)
                newpath = '%s%s/%s'%(lib, str(dircount), f)
                print('%s: mv %s %s'%(count, oldpath, newpath))
                if args.real:
                    os.rename(oldpath, newpath)
                count = 0
