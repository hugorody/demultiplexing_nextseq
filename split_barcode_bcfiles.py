#!/usr/bin/python

import sys
import time

print "Start of run:",time.ctime()

file1 = sys.argv[1]
mydir = sys.argv[2]

set_1 = set(line.strip() for line in open (file1, 'r'))

output1a = open(mydir+"bar4.txt", "w")
output1b = open(mydir+"bar5.txt", "w")
output1c = open(mydir+"bar6.txt", "w")
output1d = open(mydir+"bar7.txt", "w")
output1e = open(mydir+"bar8.txt", "w")
output1f = open(mydir+"bar9.txt", "w")


seqs = {}

for line in set_1:
	line = line.split("\t")
	a = line[0]
	b = line[1]
	seqs[a] = b


for line1 in seqs.items():
	barid = line1[0]
	barsq = line1[1]
	
	barsqsize = len(barsq)
	
	if barsqsize == 4:
		output1a.write(barid+"\t"+barsq+"\n"); #write in output1a

	if barsqsize == 5:
		output1b.write(barid+"\t"+barsq+"\n"); #write in output1b

	if barsqsize == 6:
		output1c.write(barid+"\t"+barsq+"\n"); #write in output1c

	if barsqsize == 7:
		output1d.write(barid+"\t"+barsq+"\n"); #write in output1d

	if barsqsize == 8:
		output1e.write(barid+"\t"+barsq+"\n"); #write in output1e

	if barsqsize == 9:
		output1f.write(barid+"\t"+barsq+"\n"); #write in output1f

output1a.close()
output1b.close()
output1c.close()
output1d.close()
output1e.close()
output1f.close()



individuals = set(seqs.keys())

outputind = open(mydir+"list_individuals.txt", "w")

for i in individuals:
	outputind.write(i+"\n");

outputind.close()
