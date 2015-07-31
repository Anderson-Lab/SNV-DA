import sys 


input = [row.strip("\n").split("\t") for row in open(sys.argv[1], "r")]
output = open(sys.argv[2] + "_paper.table", "w")

for row in input:
	line = ""
	line += row[0] + "&"
	line += row[3] + "&"
	line += row[1][3:] + "&"
	line += row[2][0] + "-\\(>\\)" + row[2][3] + "&"
	line += row[4] + "\\\\"
	output.write(line + "\n")	
