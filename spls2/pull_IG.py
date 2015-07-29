import csv

IgH_chr = 14
IgH_start = 106053226
IgH_end = 106518932
IgL_chr = 22
IgL_start = 22712136
IgL_end = 23248832
IgK_chr = 2
IgK_start = 106209389
IgK_end = 106725469

SNV_data = [row for row in csv.reader(open("SNVM.unsup_filtered.csv","r"))]
SNV_loci = [row for row in csv.reader(open("SNVM.unsup_loci.csv", "r"))]
output1 = csv.writer(open("SNVM.IG_data.csv", "w"))
output2 = csv.writer(open("SNVM.IG_loci.csv", "w"))

for i in range(len(SNV_data)):
	SNV = SNV_loci[i][0].split("_")
	loci = SNV[1]
	mut = SNV[2]
	chr = loci.split(":")[0][3:]
	if chr not in ["X","Y","M"]:
		chr = int(chr)
	else:
		continue
	pos = int(loci.split(":")[1])
	
	if chr == IgH_chr:
		if pos >= IgH_start and pos <= IgH_end:
			output1.writerow(SNV_data[i])
			output2.writerow(["IgH_chr" + str(chr) + ":" + str(pos) + "_" + mut, SNV_loci[i][1]])
	elif chr == IgL_chr:
        	if pos >= IgL_start and pos <= IgL_end:
                        output1.writerow(SNV_data[i])
                        output2.writerow(["IgL_chr" + str(chr) + ":" + str(pos) + "_" + mut, SNV_loci[i][1]])
	elif chr == IgK_chr:
        	if pos >= IgK_start and pos <= IgK_end:
                        output1.writerow(SNV_data[i])
                        output2.writerow(["IgK_chr" + str(chr) + ":" + str(pos) + "_" + mut, SNV_loci[i][1]])
	else:
		continue


	



