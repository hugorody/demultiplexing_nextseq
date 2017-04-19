#!/usr/bin/python
#Read filter FASTq file based on a BLAST file result
#usage: python script.py file.fastq file.blastn
 
import sys
 
filefastq = sys.argv[1]
fileblast = sys.argv[2]
 
#open files

try:
    fgen = open(filefastq)
    fbla = open(fileblast)
except IOError:
    print ("Files doesn't exist!")


reads = []     #this is a black list. Headers coming to this list will be excluded from FASTq file. Maybe change this list by tupple? I think tupple is faster.

for line in fbla:
	line = line.split("\t")
	read = line[0]
	iden = float(line[2])		#percentage of identity
	alle = int(line[3])			#alignment length
	evlu = float(line[10])		#evalue

	if iden >= 90.0 and alle >= 140 and evlu < 0.001:	#parameters the blast result should pass to be included in the reads list
		reads.append(read)

blacklist = set(reads)

#Create dictionary containing only the fastq sequences which headers are not in the reads list 
seqs={}
parsing = 0
for i in fgen:
    i = i.rstrip()
     
    if "@" in i[0]:
        scafid = i.split()[0][1:]
     
    if parsing > 0 and "@" not in i[0]:            #Second 2: if parsing is greater than 0, than this is the first line of fastq belonging to a read selected in the first step. When it finds a line with "@", then parsing will be set to 0 and loop starts again
        seqs[scafid] = seqs[scafid] + i
    else:
        parsing = 0
 
    if "@" in i[0] and scafid not in blacklist:   #First 1: if id read are not in my reads list, then parsing will be set 1
        seqs[scafid] = ''                    	   #create a dictionary with the fastq sequence correspondent to the read not in reads list
        parsing += 1


#prints the sequences from seqs dictionary
for i in seqs.items():
    header= i[0]
    seq = i[1]
    print "@"+header+"\n"+seq.replace("+", "\n+\n")
