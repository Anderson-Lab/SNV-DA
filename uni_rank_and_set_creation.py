import os
import csv
import shutil
import random
import sys

def create_sets(first_index, second_index):

	featureRows = {}
	testRows = {}
	featureValues = {}
	i = 0
	thrown_out= 0
	for x in features:
		snps = covered[x]
		if snps[first] == "Na" and snps[second] == "Na":
							thrown_out += 1
							continue			


		masked = []
		
		sumSNPs = float(0)
		amount = 0

		knownFIRST = 0
		FIRSTsnps = 0
		sumFIRST = 0
		
		knownSEC = 0
		SECsnps = 0
		sumSEC = 0
		FIRSTsnpsNA = 0
		SECsnpsNA = 0
		sumFIRST_NA = 0
		sumSEC_NA = 0
		
		for i in range(0,len(snps),1):
			
			if snps[i] != "Na":
				sumSNPs += float(snps[i])
				amount += 1
				if i < firstSize:
					knownFIRST += 1
					if float(snps[i]) > 0 and i != first and i != second:
						FIRSTsnps += 1
						sumFIRST += float(snps[i])
					sumFIRST_NA += float(snps[i]) 
				else:
					knownSEC += 1
					if float(snps[i]) > 0 and i != first and i != second:
						SECsnps += 1
						sumSEC += float(snps[i])
					sumSEC_NA += float(snps[i])
			else:
				masked.append(i)

		
		if FIRSTsnps != 0:
			avgFIRST = float("{0:.5f}".format(sumFIRST/FIRSTsnps))
		else:
			avgFIRST = 0
		
		if SECsnps != 0:
			avgSEC = float("{0:.5f}".format(sumSEC/SECsnps))
		else:
			avgSEC = 0


		if knownFIRST != 0:
			avgFIRST_NA = float("{0:.5f}".format(sumFIRST_NA/knownFIRST))
		else:
			avgFIRST_NA = 0
					
		if knownSEC != 0:
			avgSEC_NA = float("{0:.5f}".format(sumSEC_NA/knownSEC))
		else:
			avgSEC_NA = 0

		NA = (avgSEC_NA + avgFIRST_NA)/float(2)
		
		absDiff = abs(float(FIRSTsnps)-float(SECsnps))
		
		if FIRSTsnps > SECsnps:
			featureValues[x] = (absDiff, avgFIRST)
		else:
			featureValues[x] = (absDiff, avgSEC)



		row = []
		testRow = []

		if snps[first] == "Na":
			testRow.append(NA)
		else:
			testRow.append(snps[first])
		
		if snps[second] == "Na":
			testRow.append(NA)
		else:
			testRow.append(snps[second])
		
		testRows[x] = testRow

		for i in range(0,len(snps),1):
			if i != first and i != (second):
				if i in masked:
					row.append(NA)
				else:
					row.append(snps[i])
		
	
		featureRows[x] = row
	

	values = featureValues.items()
	values.sort(reverse=True, key=lambda x:(x[1]))
	
	print str(first) + "_" + str(second-firstSize)
	for w in sizes:
		if w <= len(features)-thrown_out:
			name = outputName + "_top" + str(w) + "_diff"
			
			os.chdir(name)
			topFeatures = []
			
			for a in range(w):
				topFeatures.append(values[a][0])

			trainingMatrix = csv.writer(open("training.data." + str(first) + "." + str(second-11) + ".csv", "w"))
			testMatrix = csv.writer(open("test.data." + str(first) + "." + str(second-11) + ".csv", "w"))
			featureNames = csv.writer(open("feature.names." + str(first) + "." + str(second-11) + ".csv", "w"))

			for x in topFeatures:
				trainingMatrix.writerow(featureRows[x])
				testMatrix.writerow(testRows[x])
				featureNames.writerow([x])
			os.chdir("..")
	return



SNPM = csv.reader(open(sys.argv[1],"r"))
firstSize = int(sys.argv[2])
secSize = len(SNPM.next()) - (firstSize + 2)
sizes = [int(i) for i in sys.argv[3].split(",")]
annot = sys.argv[4].split(",")
outputName = sys.argv[5]
mode = sys.argv[6]

features = set()
covered = {}
for row in SNPM:
	if sum(1 for x in annot if row[1].find(x) != -1):
		covered[row[0]] = row[2:]
		features.add(row[0])

name = outputName + "_DiffSets"
if os.path.isdir(name):
	shutil.rmtree(name)
os.makedirs(name)
os.chdir(name)


for w in sizes:
	name = outputName + "_top" + str(w) + "_diff"
	os.makedirs(name)


firstRange = [int(i) for i in range(0,firstSize,1)]
secondRange = [int(i) for i in range(firstSize,firstSize+secSize, 1)]

if mode == "every_combo":

	for first in firstRange:
		for second in secondRange:
			create_sets(first, second)

elif mode == "rand_subset":
	completed = set()
	subsetNum = int(sys.argv[7])
	for i in range(0, subsetNum, 1):
		first = random.choice(firstRange)
		second = random.choice(secondRange)
		while str(first) + "_" + str(second) in completed:
			first = random.choice(firstRange)
			second = random.choice(secondRange)
	
		create_sets(first, second)
		completed.add(str(first) + "_" + str(second))




