#!/usr/bin/python
#Read filter FASTq file based on a BLAST file result
#usage: python script.py file.fastq file.blastn
 
import sys
from Bio import SeqIO

filefastq = sys.argv[1]
fileblast = sys.argv[2]

with open(fileblast,'r') as fbla:

    blacklist = {}   #this is a black list. Headers coming to this list will be excluded from FASTq file. Maybe change this list by tupple? I think tupple is faster.
    
    for line in fbla:
    	line = line.split("\t")
    	read = line[0]
    	iden = float(line[2])		#percentage of identity
    	alle = int(line[3])		#alignment length
    	evlu = float(line[10])		#evalue
    
	if iden >= 90.0 and alle >= 140 and evlu < 0.001:	#parameters the blast result should pass to be included in the reads list
    		blacklist[read] = ''

#read fastq and check if it is not on blacklist
for record in SeqIO.parse(filefastq, "fastq"):
    if record.id not in blacklist:
        print record.format("fastq") #full sequence with header, sequence, and phred quality letters
