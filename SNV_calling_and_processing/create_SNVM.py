import os
import csv
import sys

uniqs = set()
anno = {}
annos = csv.reader(open(sys.argv[3],"r"))
for row in annos:
	anno[row[0]+ ":" + row[1] + "_" + row[3] + "->" + row[4]] = (row[6],row[5],row[8])
	uniqs.add(row[0]+ ":" + row[1] + "_" + row[3] + "->" + row[4])

covDir = sys.argv[1]
samples = []
cov = {}
for file in os.listdir(covDir):
	if file.find("snv_loci_cov") != -1:
		sample = file[:file.find(".")]
		samples.append(sample)
		print sample
		cov[sample] = {}
		for x in uniqs:
			cov[sample][x] = 0
		covs = open(covDir + file,"r")
		for row in covs:
			line = row.strip("\n").split("\t")
			cov[sample][line[4] + ":" + line[6] + "_" + line[7] + "->" + line[8]] = 1

freqsDir = sys.argv[2]
freq = {}
for file in os.listdir(freqsDir):
	if file.find(".freq") != -1:
		sample = file[:file.find(".")]
		freq[sample] = {}
		for x in uniqs:
			freq[sample][x] = 0
		freqs = open(freqsDir + file,"r")
		for row in freqs:
			line = row.strip("\n").split("\t")
			freq[sample][line[0]] = line[1]


uniqs = set()
anno = {}
annos = csv.reader(open(sys.argv[3],"r"))
for row in annos:
	anno[row[0]+ ":" + row[1] + "_" + row[3] + "->" + row[4]] = (row[6],row[5],row[8])
	uniqs.add(row[0]+ ":" + row[1] + "_" + row[3] + "->" + row[4])


snpm = csv.writer(open("SNVM.unfiltered.csv","w"))
firstRow = ["SNV", "annot"]
for sample in samples:
	firstRow.append(sample)
snpm.writerow(firstRow)

for SNP in uniqs:
	row = [str.replace(anno[SNP][0],",","&") + "_" + SNP]
	annot = anno[SNP][1]
	if anno[SNP][1] == "exonic":
		if anno[SNP][2] == "nonsynonymous SNV":
			annot = "nonsyn_exonic"
		elif anno[SNP][2] == "synonymous SNV":
			annot = "syn_exonic"
		else:
			annot = "exonic"
	row.append(annot)
	for sample in samples:
		if cov[sample][SNP] == 1:
			row.append(freq[sample][SNP])
		else:
			row.append("Na")

	snpm.writerow(row)

