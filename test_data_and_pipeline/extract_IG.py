import csv

input = csv.reader(open("SNVM.unfiltered.csv", "r"))
output = csv.writer(open("Ig_snvs.csv", "w"))
names = csv.writer(open("Ig_loci.csv", "w"))
input.next()
for row in input:
	type = row[1]
	if type == "intergenic":
		pos = row[0][row[0].find("chr"):].split("_")[0]
		sub_pos = pos.split(":")
		chr = sub_pos[0][3:]
		bp = sub_pos[1]
		if (chr == "14" and int(bp) >= 106032614 and int(bp) <= 107288051) or (chr == "22" and int(bp) >= 22380474 and int(bp) <= 23265085) or (chr == "2" and int(bp) >= 89890568 and int(bp) <= 90274235):
			nonzero = 0
			for i in range(2, len(row), 1):
				if row[i] != "Na":
					if float(row[i]) > 0:
						nonzero += 1 

			if nonzero > 4:
				output.writerow(row[2:])
				names.writerow(row[:2])




