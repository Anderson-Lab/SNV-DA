import csv
import sys

input = csv.reader(open("SNVM.unfiltered.csv", "r"))
snvs = [x.strip("\n") for x in open("intronic_snvs.txt", "r")]
print snvs
output = csv.writer(open(sys.argv[1] + "_top_snvs_report.csv", "w"))
output2 = csv.writer(open(sys.argv[1] + "_top_snvs_info.csv", "w")) 
output3 = open(sys.argv[1] + "_vep.csv", "w")
infos = {}
for row in input:
	found = 0
	found = sum([1 for i in range(len(snvs)) if row[0].find(snvs[i]) != -1])
	if found > 0:
		output2.writerow(row)
		gene = row[0].split("_")[0]
		pos = row[0].split("_")[1]
		mut = row[0].split("_")[2]
		chr = pos.split(":")[0][3:]
		bp = pos.split(":")[1]
		first = mut.split("->")[0]
		sec = mut.split("->")[1]
		output3.write(chr + " " + bp + " " + first + " " + sec + "\n")
		numDF = 0
		numR = 0
		totalDF = 0
		totalR = 0
		sumDF = 0
		sumR = 0
		for i in range(20):
			if row[i+2] != "Na":
				if i < 11:
					totalDF += 1
				else:
					totalR += 1		
		


				if float(row[i+2]) > 0:
					if i < 11:
						numDF += 1
						sumDF += float(row[i+2])
					else:
						numR += 1
						sumR += float(row[i+2])
		avgDF = 0
		avgR = 0
		if numDF > 0:
			avgDF = sumDF / float(numDF)
		if numR > 0:
			avgR = sumR / float(numR)
		infos[row[0]] = [gene, "chr" + chr + ":" + bp, mut,  numDF, numR, float(numDF) / float(totalDF), float(numR) / float(totalR), avgDF, avgR]

print infos.keys()
for x in snvs:
	output.writerow(infos[x])
