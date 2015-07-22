import csv
import math

input = csv.reader(open("SNVM.unfiltered.csv","r"))
input.next()
indy = 0
loci = {}
for row in input:
	SNVloci = row[0][:-5]
	snv = row[0].strip('"')
	#print snv	
	if SNVloci not in loci.keys():
		loci[SNVloci] = [[snv[len(snv)-4:]]+row[1:]]
	else:
		loci[SNVloci] += [[snv[len(snv)-4:]]+row[1:]]
	print indy
	indy += 1
#	if indy == 2000:
#		break

output = csv.writer(open("SNVM.fishers_rank.csv","w"))
snv_ps = []

for key in loci.keys():
	tab = [[0,0,0,0], [0,0,0,0]]
	snvs = loci[key]
	base = snvs[0][0][0]
	if base == "A":
		bindex = 0
	elif base == "C":
		bindex = 1
	elif base == "G":
		bindex = 2
	else:
		bindex = 3

	base_df = 11
        base_r = 9

	base_amounts = [1 for x in range(20)]
	for i in range(len(snvs)):
		snv = snvs[i]
		bp = snv[0][3]
		print bp
		if bp == "A":
                	aindex = 0
        	elif bp == "C":
                	aindex = 1
        	elif bp == "G":
                	aindex = 2
        	else:
                	aindex = 3
		for z in range(2, len(snv),1):
			class_index = 0
			if z > 12:
				class_index = 1

			if snv[z] == "Na":
				base_amounts[z-2] = 0
			else:
				if float(snv[z]) > 0:
					tab[class_index][aindex] += 1
					base_amounts[z-2] -= float(snv[z])
	
	
	tab[0][bindex] = sum([1 for y in range(11) if float(base_amounts[y]) > 0])
	tab[1][bindex] = sum([1 for y in range(11, 20,1) if float(base_amounts[y]) > 0])
	As = tab[0][0] + tab[1][0]
	cs = tab[0][1] + tab[1][1]
	gs = tab[0][2] + tab[1][2]
	ts = tab[0][3] + tab[1][3]
	dfs = sum(tab[0])
	rs = sum(tab[1])
	n = dfs + rs
	p = float(math.factorial(As) * math.factorial(cs) * math.factorial(gs) * math.factorial(ts) * math.factorial(dfs) * math.factorial(rs)) 
	denom = 1
	for r in range(3):
		for c in range(2):
			denom = denom * float(math.factorial(tab[c][r]))
	denom = denom * math.factorial(n)
	p = p / float(denom)
	snv_ps.append([key,p, tab, snv[1]])

snv_ps.sort(key=lambda tup: tup[1])
for x in snv_ps:
	output.writerow(x)



