#!/bin/bash/python

import subprocess

ndirs = 11
dirs = [ (str(i) + '/demux/') for i in range(ndirs) ]
base = '/fh/fast/bedford_t/zika-seq/data/usvi-library7-1d-2017-03-24/test/workspace/'
out = base + 'demux/'
demux_gzs = [ ( 'NB%02d.fasta.gz' % (i) ) for i in range(1,13) ]
demux_fastas = [ ( 'NB%02d.fasta' % (i) ) for i in range(1,13) ]

for d in dirs:
    path = base + d
    for g in demux_gzs:
	g1 = path + g
	call = 'gunzip %s' % (g1)
	print(call)
	subprocess.call(call, shell=True)
    print('\n#####\n\n')

nb01s = [ base + o + demux_fastas[0] for o in dirs ]
nb02s = [ base + o + demux_fastas[1] for o in dirs ]
nb03s = [ base + o + demux_fastas[2] for o in dirs ]
nb04s = [ base + o + demux_fastas[3] for o in dirs ]
nb05s = [ base + o + demux_fastas[4] for o in dirs ]
nb06s = [ base + o + demux_fastas[5] for o in dirs ]
nb07s = [ base + o + demux_fastas[6] for o in dirs ]
nb08s = [ base + o + demux_fastas[7] for o in dirs ]
nb09s = [ base + o + demux_fastas[8] for o in dirs ]
nb10s = [ base + o + demux_fastas[9] for o in dirs ]
nb11s = [ base + o + demux_fastas[10] for o in dirs ]
nb12s = [ base + o + demux_fastas[11] for o in dirs ]

print('nb01s: ' + " ".join(nb01s))
print('nb02s: ' + " ".join(nb02s))
print('nb03s: ' + " ".join(nb03s))
print('nb04s: ' + " ".join(nb04s))
print('nb05s: ' + " ".join(nb05s))
print('nb06s: ' + " ".join(nb06s))
print('nb07s: ' + " ".join(nb07s))
print('nb08s: ' + " ".join(nb08s))
print('nb09s: ' + " ".join(nb09s))
print('nb10s: ' + " ".join(nb10s))
print('nb11s: ' + " ".join(nb11s))
print('nb12s: ' + " ".join(nb12s))

print('\n#####\n\n')

barcodes = [ nb01s, nb02s, nb03s, nb04s, nb05s, nb06s, nb07s, nb08s, nb09s, nb10s, nb11s, nb12s ]

counter = 1
for b in barcodes:
    fname = '%sNB%02d_complete.fasta' % (out,counter)
    call = " ".join(['cat'] + b) + ' > %s' % (fname)
    f = open(fname, 'w+')
    print(fname)
    print(call)
    subprocess.call(call, shell=True, stdout=f)
    counter += 1
