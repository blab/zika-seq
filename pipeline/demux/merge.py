#!/bin/bash/python

import os, subprocess
import time

t = time.time()
count = 0
cwd = os.getcwd()
print("cwd: %s" % cwd)

for folder in os.listdir(cwd):
    path = os.path.join(cwd, folder)
    print('~~~~~ path: %s ~~~~~' % (path))
    if os.path.isdir(path):
        for f in os.listdir(path):
            call = 'mv %s .' % (path + '/' + f)
            print(call)
            subprocess.call(call, shell=True)
            count += 1

        call = 'rmdir %s' % (path)
        print(call)
        subprocess.call(call, shell=True)

print("moved %s files in %s seconds" % (count, (time.time()-t)))
