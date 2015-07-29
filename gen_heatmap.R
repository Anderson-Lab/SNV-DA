library('ggplot2')
pdf("nonsyn_snvs_heatmap.pdf", 7,5)
input = read.csv("nonsyn_snp.info2", header = F, na.strings="Na")
SNVs = input[1]
data = input[,3:ncol(input)]
samps = c("DF0", "DF1", "DF2","DF3", "DF4", "DF5","DF6", "DF7", "DF8","DF9", "DF10", "R0","R1","R2","R3","R4","R5","R6","R7","R8")
colnames(data) = samps
rownames(data) = t(SNVs)
data = matrix(as.numeric(unlist(data)), nrow=nrow(data))
colnames(data) = samps
rownames(data) = t(SNVs)
library(gplots)
heatmap.2(data, key.title="NA", key.xlab="Allele Fraction",density.info='none',na.color='black', col = colorpanel(100,"white","blue"),margins=c(4.2,10.4),trace='none')
dev.off()

pdf("intronic_heatmap.pdf", 7,5)
input = read.csv("intronic_snvs.info2", header = F, na.strings="Na")
SNVs = input[1]
data = input[,3:ncol(input)]
samps = c("DF0", "DF1", "DF2","DF3", "DF4", "DF5","DF6", "DF7", "DF8","DF9", "DF10", "R0","R1","R2","R3","R4","R5","R6","R7","R8")
colnames(data) = samps
rownames(data) = t(SNVs)
data = matrix(as.numeric(unlist(data)), nrow=nrow(data))
colnames(data) = samps
rownames(data) = t(SNVs)
library(gplots)
heatmap.2(data, key.title="NA", key.xlab="Allele Fraction",density.info='none',na.color='black', col = colorpanel(100,"white","blue"),margins=c(4.2,11.8),trace='none')
dev.off()