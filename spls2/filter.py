import csv

input = csv.reader(open("SNVM.unfiltered.csv", "r"))
output = csv.writer(open("SNVM.unsup_filtered.csv", "w"))
names = csv.writer(open("SNVM.unsup_loci.csv", "w"))
output2 = csv.writer(open("SNVM.unsup_total.csv", "w"))
input.next()
for row in input:
	type = row[1]
	class1 = 0
	class2 = 0
	for i in range(2, 2 + 11, 1):
		if row[i] != "Na":
			if float(row[i]) > 0:
				class1 += 1 

	for i in range(13, len(row), 1):
                if row[i] != "Na":
                	if float(row[i]) > 0:
                        	class2 += 1 
	if class1 + class2 > 2:
		output.writerow(row[2:])
		names.writerow(row[:2])
		output2.writerow(row)




