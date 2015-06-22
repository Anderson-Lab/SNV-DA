
library(mixOmics)
args <- commandArgs(trailingOnly = TRUE)
#args = c("exonic", "test2")


allSNPs <- data.matrix(read.csv("~/repos/SNV-DA/test_data_and_pipeline/SNVM.filtered.csv", header=T, stringsAsFactors=F,sep=","))
names <- read.csv("~/repos/SNV-DA/test_data_and_pipeline/SNV_info.csv", header=F, sep=',')
rownames(allSNPs) = names[,1]
#allSNPs = allSNPs[1:300,]


types = strsplit(args[1], ",")[[1]]

if (types[1] != "none"){
	typed = c()
	for (i in 1:length(types)){

		typed = union(rownames(allSNPs)[which(grepl(types[i],names[,2]))], typed) 
	}
}else{
	typed = rownames(allSNPs)
}


snvs = rownames(allSNPs)
num_class_1 = 11
name = args[2]
testAmounts = seq(5, 500, 15)
accuracy = c()

for (numFeatures in testAmounts){
	print(numFeatures)
	correct <- 0
	missed <- 0
	iter = 0
	for (skip1 in 1:num_class_1){
		for (skip2 in (num_class_1+1):ncol(allSNPs)){
				iter = iter + 1
				print(iter)
				
				test.data = na.omit(allSNPs[,c(skip1, skip2)])
				keep = rownames(test.data)
				training.data = allSNPs[,c(-skip1,-skip2)][which(rownames(allSNPs) %in% keep & rownames(allSNPs) %in% typed ),]
				head(training.data)
				splsda.model <- splsda(t(training.data), factor(t(c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2))), ncomp=1, keepX=numFeatures)
				training.data = training.data[which(rownames(training.data) %in% colnames(splsda.model$X)),]
	# 			print("~")
	# 			print(nrow(training.data))
	# 			print("!")
				test.data = test.data[which(rownames(test.data) %in% colnames(splsda.model$X) & rownames(test.data) %in% keep  & rownames(training.data) %in% typed),]
	# 			print(nrow(test.data))
	# 			print("~")
				prediction <- predict(splsda.model, t(test.data))
				pred.classes <- prediction$class$centroids.dist;
				pred.class.one <- pred.classes[1];
				pred.class.two <- pred.classes[2];
				
				if (pred.class.one == 1){
					correct <- correct + 1;
				}else{
					missed <- missed + 1;
				}
				
				if (pred.class.two == 2){
					correct <- correct + 1;
				}else{
					missed <- missed + 1;
				}
				accuracy = c(accuracy, correct / (correct + missed))
			
	}
	}
}  
write(t(cbind(c(testAmounts), c(accuracy)), ncol=2,file=paste(name, "_cv_accs.txt", sep=""), sep="\t")
