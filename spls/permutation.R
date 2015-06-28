
library(mixOmics)
library(pROC)
args <- commandArgs(trailingOnly = TRUE)
name = args[4]
##args = c("none", "test2", 5, 15, 5)


allSNPs <- data.matrix(read.csv(args[1], header=F, stringsAsFactors=F,sep=","))
names <- read.csv(args[2], header=F, sep=',')
rownames(allSNPs) = names[,1]
#allSNPs = allSNPs[1:300,]


types = strsplit(args[3], ",")[[1]]

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



optM = as.numeric(args[5])
numIter = args[6]
AUCs = c()
for (run in 1:numIter){
####Final CV using the average of the bestM's
correct <- 0
missed <- 0
classes <- c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)
indices <- sample(c(seq(1,20,1)))
first_class <- indices[1:num_class_1]
second_class <- indices[num_class_1+1:length(indices)]
second_class = second_class[!is.na(second_class)]

iter = 0
resp <- c()
pred <- c()
class1 <- c()
class2 <- c()
print(first_class)
print(second_class)
for (skip1 in first_class){
  for (skip2 in second_class){
    iter = iter + 1
    
    test.data = na.omit(allSNPs[,c(skip1, skip2)])
    keep = rownames(test.data)
    training.data = allSNPs[,c(-skip1,-skip2)][which(rownames(allSNPs) %in% keep & rownames(allSNPs) %in% typed ),]
    head(training.data)
    splsda.model <- splsda(t(training.data), factor(t(classes)), ncomp=1, keepX=optM)
    training.data = training.data[which(rownames(training.data) %in% colnames(splsda.model$X)),]

    test.data = test.data[which(rownames(test.data) %in% colnames(splsda.model$X)),]

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
    #Computing roc values
    unrounded.prediction.1 <- prediction$variates[1]
    unrounded.prediction.2 <- prediction$variates[2]
    
    class.1.centroid <- prediction$centroid[1]
    class.2.centroid <- prediction$centroid[2]
    
    inv <- (class.1.centroid > class.2.centroid)
    
    if(inv){
      temp <- class.1.centroid
      class.1.centroid <- class.2.centroid
      class.2.centroid <- temp
    }
    
    class.2.centroid <- class.2.centroid - class.1.centroid
    unrounded.prediction.1 <- unrounded.prediction.1 - class.1.centroid
    unrounded.prediction.2 <- unrounded.prediction.2 - class.1.centroid
    class.1.centroid <- 0
    
    unrounded.prediction.1 <- unrounded.prediction.1 / class.2.centroid
    unrounded.prediction.2 <- unrounded.prediction.2 / class.2.centroid
    class.2.centroid <- 1
    
    if(inv){
      unrounded.prediction.1 <- 1 - unrounded.prediction.1
      unrounded.prediction.2 <- 1 - unrounded.prediction.2
    }
    
    
    class1 = c(class1, unrounded.prediction.1)
    class2 = c(class2, unrounded.prediction.2)
  }	
}

print(indices)


auc_ci = as.numeric(ci.auc(roc(controls=class1, cases=class2)))
AUCs = c(AUCs, auc_ci[2])
print(run)
}
write(t(AUCs), ncol=1, file=paste(name, "_permutation.txt", sep=""), sep="\t")
