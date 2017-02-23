import argparse
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory', default=None, help="Global path to directory of native barcode base-called .fast5's")

def preprocess(path):
    '''
    This function converts .fast5 files in path/NB*/ into input for depth_coverage.R
    Steps are:
      1. poretools fasta --type 2D <path/to/base/called/reads/> > <name.fasta>
      2. bwa mem -x ont2d <indexed_reference.fasta> <name.fasta> | samtools view -bS - | samtools sort -o <name.sorted.bam> -
      3. samtools depth <name.sorted.bam> > <name.coverage>
      3a.head <name.coverage> # This finds the name of the 'chromosome'; there may be >1
      4. awk '$1 == "<chromosomename>" {print $0}' <name.coverage> > chr1.coverage
      5. Repeat for all libraries
    '''
    for subdir in os.listdir(path):

        print("Processing library "+subdir)
        prdir = path+subdir+"/processed/"
        if not os.path.exists(prdir):
            print("Making directory "+prdir)
            os.makedirs(prdir)
        fName = prdir+subdir+".fasta"
        print("Building .fasta from .fast5's in "+subdir)
        fCall = "poretools fasta --type 2D "+path+subdir+" > "+fName
        print("$ "+fCall) #Test
        subprocess.call(fCall, shell=True)

        print("Done building "+subdir+".fasta, beginning bwa mem")
        bCall = "bwa mem -x ont2d refs/ZikaReferenceGenome.fasta "+fName+" | samtools view -bS - | samtools sort -o "+prdir+subdir+".sorted.bam -"
        subprocess.call(bCall, shell=True)
        print("$ "+bCall) #Test

        print("Done with bwa mem, beginning samtools depth")
        dCall = "samtools depth "+prdir+subdir+".sorted.bam > "+prdir+subdir+".coverage"
        # tCall = "head "+prdir+subdir+".coverage"
        subprocess.call(dCall, shell=True)
        print("$ "+dCall) #Test

        print("Done with samtools depth, beginning awk")
        aCall = "awk '$1 == \"NC_012532.1\" {print $0}' "+prdir+subdir+".coverage > "+prdir+"chr1.coverage"
        subprocess.call(aCall, shell=True)
        print("$ "+aCall) #Test

if __name__=="__main__":

    params = parser.parse_args()

    p = params.directory
    preprocess(p)
