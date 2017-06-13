#!/usr/bin/python

import sys
from Bio import SeqIO
import glob, os

file1 = sys.argv[1] #list containing individuals ids separated by breaklines
dir1 = sys.argv[2] #input directory, where the fq files
dir2 = sys.argv[3] #output directory, to write the fastq files

#create a list of individuals based on file1
listindividuals = []
with open(file1) as set1:
    for i in set1:
        i = i.rstrip()
        listindividuals.append(i)

#for each in list of individuals
for i in listindividuals:

    searchfile = i+".fq"
    outputname = str(dir2)+str(i)+".fastq" #outputname is set to write in the output directory 2
    output = open(outputname,'w')

    lisftoffiles = [] #list of files that match to individual
    
    os.chdir(dir1) #list all files in the input directory 1
    for file in glob.glob("*.fq"): #and filter by *.fq
        if searchfile in file: #if search file, which is with .fq, is on any of the files from input dir 1
            lisftoffiles.append(file) #append the name of the file in list of files
    
    for j in lisftoffiles: #for each of the files in list of files that match to the individual
        for record in SeqIO.parse(j, "fastq"): #read each read in fastq format
            output.write(record.format("fastq")); #write to the output named individual.fastq
    output.close() #close the output
