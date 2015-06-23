
library(mixOmics)
library(pROC)
args <- commandArgs(trailingOnly = TRUE)
name = args[2]
##args = c("none", "test2", 5, 15, 5)


allSNPs <- data.matrix(read.csv("~/repos/SNV-DA/test_data_and_pipeline/Ig_snvs.csv", header=F, stringsAsFactors=F,sep=","))
names <- read.csv("~/repos/SNV-DA/test_data_and_pipeline/Ig_loci.csv", header=F, sep=',')
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

print("Number of SNVs")
print(length(typed))
snvs = rownames(allSNPs)
num_class_1 = 11

testAmounts = seq(as.numeric(args[3]), as.numeric(args[4]), as.numeric(args[5]))
accuracy = c()
all_roc_data = matrix(,nrow=100,ncol=2*length(testAmounts))
header = c()
for(x in testAmounts){
  header = c(header, x, "")
}


AUC = c()

print(testAmounts)
###Determining optimal number of features by iteratively finding optimal bestM for each training set
bestM = c()
classes = c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2)
iter = 0
for (skip1 in 1:num_class_1){
	for (skip2 in (num_class_1+1):ncol(allSNPs)){
    iter = iter + 1
    print("Iteration: ")
    print(iter)
	  test.data = na.omit(allSNPs[,c(skip1, skip2)])
	  keep = rownames(test.data)
	  training.data = allSNPs[,c(-skip1,-skip2)][which(rownames(allSNPs) %in% keep & rownames(allSNPs) %in% typed ),]
    accuracy = c()
    for (numFeatures in testAmounts){
      print("M: ")
      print(numFeatures)
      correct = 0
      missed = 0
      for (sub_skip in 1:ncol(training.data)){
        sub.test.class = classes[sub_skip]
        sub.classes = classes[-sub_skip]
        sub.test.data = na.omit(training.data[,sub_skip])
        keep = attributes(sub.test.data)$names
        sub.training.data = training.data[,-sub_skip][which(rownames(training.data) %in% keep),]
        splsda.model <- splsda(t(sub.training.data), factor(t(sub.classes)), ncomp=1, keepX=numFeatures)
        #sub.training.data = sub.training.data[which(rownames(sub.training.data) %in% colnames(splsda.model$X)),]
				sub.test.data = sub.test.data[which(attributes(sub.test.data)$names %in% colnames(splsda.model$X))]
				prediction <- predict(splsda.model, t(sub.test.data))
				pred.classes <- prediction$class$centroids.dist;
				pred.class.one <- pred.classes[1];

				if (pred.class.one == sub.test.class){
					correct <- correct + 1;
				}else{
					missed <- missed + 1;
				}
      }
    accuracy = c(accuracy, correct / (missed + correct)) 
    }

  bestM = c(bestM, round(mean(testAmounts[which(accuracy == max(accuracy))])))
  print(bestM)
  }
}
write(t(c(bestM)), ncol=1, file=paste(name, "_iteration_optimals.txt", sep=""), sep="\t")




####Final CV using the average of the bestM's
correct <- 0
missed <- 0

iter = 0
resp <- c()
pred <- c()
class1 <- c()
class2 <- c()
optM = round(mean(bestM))
for (skip1 in 1:num_class_1){
  for (skip2 in (num_class_1+1):ncol(allSNPs)){
    iter = iter + 1
    print(iter)
    
    test.data = na.omit(allSNPs[,c(skip1, skip2)])
    keep = rownames(test.data)
    training.data = allSNPs[,c(-skip1,-skip2)][which(rownames(allSNPs) %in% keep & rownames(allSNPs) %in% typed ),]
    head(training.data)
    splsda.model <- splsda(t(training.data), factor(t(classes)), ncomp=1, keepX=optM)
    training.data = training.data[which(rownames(training.data) %in% colnames(splsda.model$X)),]

    test.data = test.data[which(rownames(test.data) %in% colnames(splsda.model$X) & rownames(test.data) %in% keep  & rownames(training.data) %in% typed),]

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

    #Computing ROC values
    unrounded.prediction.1 <- prediction$variates[1]

    unrounded.prediction.2 <- prediction$variates[2]

    
    class1 = c(class1, unrounded.prediction.1)
    class2 = c(class2, unrounded.prediction.2)
  }	
}



print("Optimal M = ")
print(optM)

auc_ci = as.numeric(ci.auc(roc(controls=class1, cases=class2)))
print("AUC CI = ")
print(auc_ci)

accuracy =  correct / (correct + missed)
print(accuracy)		
stats = c(optM, accuracy, auc_ci) 
write(stats, ncol=5, file=paste(name, "_performance_stats.txt", sep=""), sep="\t")

##Determine most relevant biomarkers with optimal M
allData = allSNPs[which(rownames(allSNPs) %in% typed ),]
total_classes = c(rep(0, num_class_1), rep(1, ncol(allData) - num_class_1))
splsda.model <- splsda(t(allData), factor(t(total_classes)), ncomp=1, keepX=optM)
SNVs = splsda.model$loadings$X
all_snv_names = rownames(SNVs)
snvs = matrix(, nrow=optM, ncol=3)
index = 1
for(i in 1:length(SNVs)){
  if (abs(SNVs[i]) > 0){
    snvs[index,1:2] = c(all_snv_names[i], abs(SNVs[i]))
    index = index + 1
  }
}
snvs[,3] = as.character(names[,2][which(names[,1] %in% snvs[,1] )])


write(t(snvs[order(snvs[,2], decreasing = T),]), ncol=3, file=paste(name, "_top_SNVs.txt", sep=""), sep="\t")




