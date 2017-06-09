#!/usr/bin/python
#usage: python script.py file1.fq file2.fq ...[all fastq files]... output.fq
#the last element of sys.argv is the output file

import sys
from Bio import SeqIO

outputname = sys.argv[len(sys.argv)-1]
output = open(outputname,'w')

for i in sys.argv[1:-1]:
        for record in SeqIO.parse(i, "fastq"):
                output.write(record.format("fastq"));

output.close()
