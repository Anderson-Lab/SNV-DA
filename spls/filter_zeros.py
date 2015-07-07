import csv

input = csv.reader(open("SNVM.unfiltered.csv", "r"))
output = csv.writer(open("SNVM.total.csv","w"))
loci = csv.writer(open("SNVM.total_loci.csv", "w"))
input.next()

for row in input:
	
	for i in range(2,len(row),1):
		if row[i] != "Na":
			if float(row[i]) > 0:
				output.writerow(row[2:])
				loci.writerow(row[:2])
				break
