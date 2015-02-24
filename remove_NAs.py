import csv
import sys

snvm = csv.reader(open(sys.argv[1], "r"))
newSnvm = csv.writer(open("SNVM.filtered.csv","w"))

newSnvm.writerow(snvm.next())

firstSize = int(sys.argv[2])
firstNa = int(sys.argv[3])
secNa = int(sys.argv[4])

for row in snvm:
	if sum(1 for x in row[2:firstSize+2] if x == "Na") <= firstNa and sum(1 for x in row[firstSize+2:len(row)] if x == "Na") <= secNa:
		newSnvm.writerow(row)
