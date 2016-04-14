import sys
import os

variantDir = sys.argv[1]
coverage = sys.argv[2]

uniqueSNPs = set()

os.system("mkdir freqs")
for file in os.listdir(variantDir):
	if file.find("variant") != -1:
		
		input = open(variantDir + "/" + file, "r")
		print variantDir + "/" + file
		freqFile = open("freqs/" + file[:file.find(".variant")] + ".freq", "w")
		for row in input:
			line = row.strip("\n").split("\t")
			if int(line[3].split(",")[0]) >= int(coverage):
				uniqueSNPs.add((line[0],int(line[2]), line[4], line[5]))
				freqFile.write(line[0]+ ":" + line[2] + "_" + line[4] + "->" + line[5] + "\t" + line[6] + "\n")

print "Number of unique SNPs with coverage greater >= 10:"
print len(uniqueSNPs)

uniq = open("high_cov_unique_SNVs.bed", "w")
anno = open("high_cov_unique_SNVs.anno_input", "w")
for x in uniqueSNPs:
	uniq.write(x[0] + "\t" + str(x[1]-1) + "\t" + str(x[1]) + "\t" + x[2] + "\t" + x[3] + "\n")
	anno.write(x[0] + "\t" + str(x[1]) + "\t" + str(x[1]) + "\t" + x[2] + "\t" + x[3] + "\n")
