import os
import sys

cov = sys.argv[1]
inName = sys.argv[2]
outName = sys.argv[3]

snvs = set()
for x in open("main.snvs", "r"):
	snvs.add(x.split("_")[1])

input = open(inName , "r")
out = open(outName, "w")
for row in input:
	split = row.split("\t") 
	numReads = int(split[2].split(",")[0])
	if numReads >= int(cov) and split[0]+":" + split[1] in snvs:
		out.write(row)
