
output = open("SNVM.just_loci.csv", "w")
for row in open("SNVM.unsup_loci.csv", "r"):
	line = row.strip("\n").split(",")
	line = line[0:len(line)-1]
	out = ""
	if len(line) > 1:
		for x in line:
			out += x + ","
		out = out[0:len(out)-1]
	else:
		out = line[0]
	output.write(out + "\n")




