#!/bin/bash/python

import os, subprocess

cwd = os.getcwd()

for folder in os.listdir(cwd):
    path = cwd + folder
    if os.path.isdir(path):
        for f in os.listdir(path):
            call = 'mv %s .'
            print(call)
            # subprocess.call(call, shell=True)
