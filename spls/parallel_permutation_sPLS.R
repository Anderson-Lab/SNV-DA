library(foreach)
library(doParallel)
library(mixOmics)
library(pROC)
args <- commandArgs(trailingOnly = TRUE)
name = args[4]
##args = c("none", "test2", 5, 15, 5)
cl<-makeCluster(12)
registerDoParallel(cl)

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

compare_auc = as.numeric(args[7])

optM = as.numeric(args[5])
numIter = as.numeric(args[6])
AUCs = c()
run = 0
#ls <- foreach (iter=1:numIter, .packages=c('mixOmics','pROC')) %dopar% {
for(iter in 1:numIter){
 write(iter, file=paste(name, ".log", sep=""), append=T)
####Final CV using the average of the bestM's
correct <- 0
missed <- 0
classes <- c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)
indices <- sample(c(seq(1,20,1)))
first_class <- indices[1:num_class_1]
second_class <- indices[num_class_1+1:length(indices)]
second_class = second_class[!is.na(second_class)]

new_allSNPs = allSNPs[which(length(intersect(which(allSNPs > 0), first_class)) > 2 | length(intersect(which(allSNPs > 0),second_class)) > 2),]

print(indices)
print(head(new_allSNPs))



run = run + 1
iter = 0
resp <- c()
pred <- c()
class1 <- c()
class2 <- c()


for (skip1 in first_class){
  for (skip2 in second_class){
    
    test.data = na.omit(new_allSNPs[,c(skip1, skip2)])
    keep = rownames(test.data)
    training.data = new_allSNPs[,c(-skip1,-skip2)][which(rownames(new_allSNPs) %in% keep & rownames(new_allSNPs) %in% typed ),]
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



auc_ci = as.numeric(ci.auc(roc(controls=class1, cases=class2)))
auc_ci[2]
}
aucs = unlist(ls)
write(unlist(ls), ncol=1, file=paste(name, "_permutation.txt", sep=""), sep="\t")
p_val = length(which(aucs >= compare_auc))/length(aucs)
print(p_val)
write(p_val, ncol=1, file=paste(name,"_p_val.txt", sep=""),sep="\t")
stopCluster(cl)
