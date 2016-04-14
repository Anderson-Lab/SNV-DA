import os
import sys

cov = sys.argv[1]
inName = sys.argv[2]
outName = sys.argv[3]


input = open(inName , "r")
out = open(outName, "w")
for row in input:
	split = row.split("\t") 
	numReads = int(split[2].split(",")[0])
	if numReads >= int(cov):
		out.write(row)
