import csv

snvs = [x.strip("\n") for x in open("SNVM.just_loci.csv", "r")]
snvInfo = {}
total = csv.reader(open("SNVM.unfiltered.csv", "r"))

for row in total:
	snvInfo[row[0]] = row[1:]

print snvInfo.keys()
output = csv.writer(open("SNVM.unsup_total.csv", "w"))
for x in snvs:
	output.writerow([x.strip("\"")] + snvInfo[x]) 
