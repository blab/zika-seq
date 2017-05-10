#!/usr/bin/env python
from Bio import SeqIO
import sys
import vcf
import subprocess
from collections import defaultdict
import os.path
import operator
from vcftagprimersites import read_bed_file

#MASKED_POSITIONS = [2282]
MASKED_POSITIONS = []

reference = sys.argv[1]
vcffile = sys.argv[2]
bamfile = sys.argv[3]

DEPTH_THRESHOLD = 20

def collect_depths(bamfile):
	if not os.path.exists(bamfile):
		raise SystemExit("bamfile %s doesn't exist" % (bamfile,))

	sys.stderr.write(bamfile)

	p = subprocess.Popen(['samtools', 'depth', bamfile],
                             stdout=subprocess.PIPE)
	out, err = p.communicate()
	depths = defaultdict(dict)
	for ln in out.split(b"\n"):
		if ln:
			contig, pos, depth = ln.split(b"\t")
			depths[contig][int(pos)] = int(depth)
	return depths

depths = collect_depths(bamfile)

def report(r, status, allele):
	idfile = os.path.basename(vcffile).split(".")[0]
	sys.stderr.write("%s\t%s\tstatus\t%s" % (idfile, r.POS, status))
	sys.stderr.write("%s\t%s\tdepth\t%s" % (idfile, r.POS, record.INFO.get('TotalReads', ['n/a'])))
	sys.stderr.write("%s\t%s\tfrac\t%s" % (idfile, r.POS, record.INFO.get('BaseCalledFraction', ['n/a'])))
	sys.stderr.write("%s\t%s\tallele\t%s" % (idfile, r.POS, allele))
	sys.stderr.write("%s\t%s\tref\t%s" % (idfile, r.POS, record.REF))

cons = ''

seq = list(SeqIO.parse(open(sys.argv[1]), "fasta"))[0]
cons = list(seq.seq)

for n, c in enumerate(cons):
	try:
		depth = depths[seq.id][n+1]
	except:
		depth = 0

	if depth < DEPTH_THRESHOLD:
		cons[n] = 'N'

for mask in MASKED_POSITIONS:
	cons[mask-1] = 'N'

sett = set()
vcf_reader = vcf.Reader(open(vcffile, 'r'))
for record in vcf_reader:
	if record.ALT[0] != '.':
		# variant call

		if record.POS in MASKED_POSITIONS:
			report(record, "masked_manual", "n")
			continue

		if 'PRIMER' in record.INFO:
			report(record, "primer_binding_site", "n")
			cons[record.POS-1] = 'N'
			continue

		support = float(record.INFO['SupportFraction'][0])
		total_reads = int(record.INFO['TotalReads'][0])
		qual = record.QUAL
		REF = record.REF
		ALT = str(record.ALT[0])

		if len(REF) > len(ALT):
			sys.stderr.write("deletion")
			continue

		if len(ALT) > len(REF):
			sys.stderr.write("insertion")
			continue

#		if support >= 0.75 and total_reads > 30:
		if qual >= 200 and total_reads >= 20 and support > 0.75:
			sys.stderr.write(REF, ALT)

			report(record, "variant", ALT)
			sett.add(record.POS)
			if len(REF) > len(ALT):
				sys.stderr.write("deletion")
				continue

			if len(ALT) > len(REF):
				sys.stderr.write("insertion")
				continue

			#for n, b in enumerate(REF):
			#	try:
			#		cons[record.POS-1+n] = ALT[n]
			#	except IndexError:
			#		cons[record.POS-1+n] = ''
			#else:
			cons[record.POS-1] = str(ALT)
		else:
			report(record, "low_qual_variant", "n")
			cons[record.POS-1] = 'N'
			continue

#sys.stderr.write(str(sett))

print(">%s" % (sys.argv[3]))
print("".join(cons))
