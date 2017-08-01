#!/bin/bash/python

import os, subprocess

cwd = os.getcwd()
print("cwd: %s" % cwd)

for folder in os.listdir(cwd):
    path = cwd + folder
    print('~~~~~ path: %s ~~~~~' % (path))
    if os.path.isdir(path):
        for f in os.listdir(path):
            call = 'mv %s .' % (path + f)
            print(call)
            # subprocess.call(call, shell=True)
