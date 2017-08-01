#!/bin/bash/python

import os,subprocess

c = os.getcwd()
call = 'mkdir demux'
print(call)
subprocess.call(call, shell=True)

for f in os.listdir(c):
    if not f.endswith('.fast5'):
        call = 'mv %s %s' % (os.path.join(c,f),os.path.join(c,'demux',path))
        print(call)
        # subprocess.call(call,shell=True)
