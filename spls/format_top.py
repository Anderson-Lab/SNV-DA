import csv
import sys

top_file = csv.reader(open(sys.argv[1],"r"))
dbSNP = open(sys.argv[3],"r")
output = csv.writer(open(sys.argv[2] + "_for_boxplot.csv","w"))
dbSNP.next()
order = []
ids = {}
for row in dbSNP:
	splits = row.strip("\n").split("\t")
	order.append(splits[0] + "_" + splits[1])
	ids[splits[0] + "_" + splits[1]] = splits[2]

top = [x for x in top_file]

for ident in order:
	for snv in top:
		if snv[0].find(ident) != -1:
			x = snv
			break
	id = "NA"
        for z in ids.keys():
        	if x[0].find(z) != -1:
                	id = ids[z]
                	break
	print x
	for i in range(2,22, 1):
		if i < 13:
			cond = "Disease Free"
		else:
			cond = "Relapse"
		if x[i] != "Na":
			print "!"
			output.writerow([id,x[i],cond])
