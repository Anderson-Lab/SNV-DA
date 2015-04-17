import os

bibsFile = open("bibs.txt", "r")
newBibs = open("newBibs.txt", "w")

bibs = {}
lines = []
for row in bibsFile:
	lines.append(row)


first = True
for i in range(0,len(lines),1):
	row = lines[i]
	if row.find("@article") != -1:
		if first != True:
			bibs[name] = (author, title, journal, year, volume, pages)
		first = False
		name = row[row.find("{")+1:row.find(",")]
	elif row.find("author") != -1:
		author = row[row.find("{")+1:row.find("}")]
	elif row.find("title") != -1:
		title = row[row.find("{")+1:row.find("}")]
	elif row.find("journal") != -1:
		journal = row[row.find("{")+1:row.find("}")]
	elif row.find("year") != -1:
		year = row[row.find("{")+1:row.find("}")]
	elif row.find("volume") != -1:
		volume = row[row.find("{")+1:row.find("}")]
	elif row.find("pages") != -1:
		pages = row[row.find("{")+1:row.find("}")]
	else:
		nothing = 1

	i += 1


for x in sorted(bibs.keys()):
	line = "\\bibitem["
	info = bibs[x]

	line += info[0][:info[0].find(",")] + " \\textit{et~al}., " + info[3] + "]{" + x + "}\n"
	
	author = info[0]
	author = author.replace(" and", "")
	newAuthor = "" + author[0]
	for c in range(1, len(author),1):
		if author[c] == " " and author[c-1] != ",":
			newAuthor += "., "
		elif author[c] == ":" :
			newAuthor += "."
		else:
			newAuthor += author[c]

	title = list(info[1])
	if title[len(title)-1] == ".":
		title[len(title)-1] = ","
	if title[len(title)-1] != ",":
		title.append(",")
	line += newAuthor + ". (" + info[3] + ") " + "".join(title) + " {\\it " + info[2] + "}, {\\bf " + info[4][:info[4].find("(")] + "}" + info[4][info[4].find("("):] + ", " + info[5] + ".\n\n"
	newBibs.write(line)




