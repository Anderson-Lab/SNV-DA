import csv
input = csv.reader(open("spls2/SNVM.unfiltered.csv", "r"))
output = open("SNV.loci", "w")

input.next()
for x in input:
        print x
	loci = x[0].split("_")[1]
        if loci[0:3] == "chr":
		chr = loci.split(":")[0]
        	
        	bp = int(loci.split(":")[1])
        else:
		loci = x[0].split("_")[2]
		chr = loci.split(":")[0]
                bp = int(loci.split(":")[1])
	output.write(chr + "\t" + str(bp) + "\t" + str(bp+1) + "\n")
