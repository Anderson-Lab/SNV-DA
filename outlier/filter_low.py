input = open("95.variants.bed.txt", "r")
output = open("Rout.loci", "w")

for row in input:
	line = row.strip("\n").split("\t")
	cov = int(line[3].split(",")[0])
	if cov > 9:
		output.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\n")	
