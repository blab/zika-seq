import shutil
import os

workspace = '/fh/fast/bedford_t/zika-seq/data/usvi-library8-1d-2017-03-31/raw_reads'
new_dir = '/fh/bpotter/test_fast5s'
count=0
for n in range(600,610):
    for f in os.listdir("%s/%s"%(workspace,n)):
        old = "%s/%s/%s" % (workspace,n,f)
        new = "%s/%s" % (new_dir,f)
        shutil.copyfile(old,new)
        count += 1

print("Moved %s files to %s." % (count,new_dir))
