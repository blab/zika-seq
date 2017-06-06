import os
from shutil import copyfile

pass_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/pass/'
fail_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/fail/'
new_dir = '/fh/fast/bedford_t/zika-seq/data/usvi-library1-2016-12-10/test/'

for barcode in os.listdir(pass_dir):
    bcd = pass_dir + barcode + '/'
    if os.path.isdir(bcd):
        for f in bcd:
            if f[:-6] == '.fast5':
                old = bcd + f
                new = new_dir + f
                copyfile(old, new)

for f in os.listdir(fail_dir):
    if f[:-6] == '.fast5':
        old = fail_dir + f
        new = new_dir + f
        copyfile(old, new)
