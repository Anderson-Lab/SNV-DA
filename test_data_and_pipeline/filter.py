import csv
import math

input = csv.reader(open("SNVM.filtered.csv", "r"))
output = csv.writer(open("SNVM.filtered2.csv", "w"))

numSamps = len(input.next())
for row in input:
	nonzero = 0
	for i in range(len(row)):
		if row[i] != "Na":
			if row[i] > 0:
				nonzero += 1


	if nonzero >= math.floor(numSamps / 5):
		output.writerow(row)




