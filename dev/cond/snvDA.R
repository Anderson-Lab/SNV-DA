suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(gplots))

option_list <- list(


 make_option(c("-M","--matrix"),type="character", default=NULL,
 help = "SNVM input (csv)"),

 make_option(c("-Z", "--numNonZeros"), type="integer", default=3,
 help="Limit SNVs to those that have at least X number of samples with allele fractions greater than zero (default = 3). Decreases runtime. Suggested value is a 1/3 the number of samples in the smaller group."),

 make_option(c("-F","--nameClass1"),type="character", default="Class_1",
 						help = "Name of class 1 (default = Class_1)."),
 
 make_option(c("-H","--nameClass2"),type="character", default="Class 2",
 						help = "Name of class 2 (default = Class_2)."),
 
 make_option(c("-D","--sizeClass1"),type="integer", default=NULL,
 help = "The number of samples of class 1. (Used to find the number of samples of both classes)."),

 make_option(c("-G", "--include"), type="character", default=NULL,
 help="SNVs with this substring are included, e.g. \"exonic\" for nonsynonymous and synonymous exonic SNVs"),

 make_option(c("-X", "--exclude"), type="character", default=NULL,
 help="SNVs with this substring are excluded, e.g. ncRNA to remove intronic_ncRNA"),

 make_option(c("-A","--name"),type="character", default="SNV-DA_run",
 help = "One word name of the analysis used in output file names and figures."),

 make_option(c("-L", "--minK"), type="integer", default=NULL,
 help="The minimum value of K tested."),

 make_option(c("-B", "--maxK"), type="integer",default=NULL,
 help="The maximum value of K tested."),

 make_option(c("-E", "--everyKth"), type="integer",default=NULL,
 help="Specifies every Nth K tested."),

 make_option(c("-R", "--range"), type="character",default=NULL,
 help="A comma seperated list of K's to be added to the Ks tested, e.g. 5,10,15,20,25,50,75,100,250,300."),

 make_option(c("-T", "--threads"), type="integer", default=1,
 help="Number of threads for cross-validations (default = 1)."),
 
 make_option(c("-CV", "--numCV"), type="integer", default=NULL,
 help="Number of cross-validations. If flagged, SNV-DA will run X number of random cross-validations without replacement."),

 make_option(c("-lo", "--numOut"), type="integer", default=1,
 help="Number of samples removed from each group during cross-validation. If --unstratified is flagged, it is the number of total samples removed during each iteration of cross-validations."),

  make_option(c("-U", "--unstratified"), action="store_true", default=F,
 help="Test samples during cross-validations will not be stratified by group. That is test samples will be from a combination of all samples."),

 make_option(c("-findOptK", "--findOptimalK"), action="store_true", default=F,
 help="If flagged, SNV-DA will run nest cross-validations to determine optimal K."),

  make_option(c("-eval", "--evaluatePerformance"), action="store_true", default=F,
 help="If flagged, SNV-DA will evaluate performance of model the value of K found from running --findOptimalK or user supplied --optimalK."),

 make_option(c("-optK", "--optimalK"), type="integer", default=500,
 help="Used when --evaluatePerformance or --permuationTests is flagged. Specifies the value of K used for cross-validation (default = 500)."),

 make_option(c("-findSNVs","--findTopKSNVs"), action="store_true", default=F, 
 help="To find top selected SNVs, trains sPLS-DA on all samples with the value of K found from running --findOptimalK or user supplied --optimalK."), 

 make_option(c("-perm", "--permutationTests"), action="store_true", default=F,
 help="If flagged, SNV-DA will run permutation tests by randomly permuting sample classes then running cross-validations. Tests are run --permIter number of times. Performances are then compared to the value of AUC found by --evaluatePerformance of user supplied --AreaUnderCurve."),

make_option(c("-permIter", "--permutationIterations"), type="integer", default=NULL,
 help="Used when --permuationTests is flagged. Specifies the number of permutation tests to be run"),
 
 make_option(c("-auc", "--AreaUnderCurve"), type="double", default=NULL,
 help="Performance of true model. Compared to permuted models when --permutationTests are run. If unflagged, value is determined by --evaluatePerformance"),

 make_option(c("-f", "--produceFigures"), action="store_true", default=F,
 help="If flagged, allele fraction box plots, heatmaps, and kernal density figures are produced.")

 )

args <- parse_args(OptionParser(option_list=option_list))
marks = c("A", "B","C","D","E","F","G","H")

#Argument handling
name = args$name

##Check to see in SNVM was input
if(is.null(args$matrix)){
	
stop("SNV matrix not provided")	
}
write(paste("Running analysis: ", name, sep=" "), file=paste(name, ".log", sep=""), append=T)
write(paste("Input SNVM:", args$matrix, sep=" "), file=paste(name, ".log", sep=""), append=T)

##Output to the log the required number of nonzero values allows
write(paste("Removing SNVs that do not have at least ", args$numNonZeros," samples with allele fraction values >0", sep=""), file=paste(name, ".log", sep=""), append=T)



##Read in SNVM and set the names of the loci. 
allSNPs <- read.csv(args$matrix, header=T, stringsAsFactors=F,sep=",")
names <- allSNPs[,c(1,2)]
allSNPs <- data.matrix(allSNPs[,c(-1,-2)])
rownames(allSNPs) = names[,1]
totalNum = nrow(allSNPs)

##Determine rows that have X non-zero feature values
if(args$numNonZeros > 1){
	indices = c()
	for(i in 1:nrow(allSNPs)){
		values = unlist(na.omit(allSNPs[i,]))
		numNZ = length(which(values > 0))
		if(numNZ >= args$numNonZeros){
			indices = c(indices, i)
		}
	}
}else{indices = rownames(allSNPs)}


##Pull out those SNVs
names = names[indices,]
allSNPs = allSNPs[indices,]
filteredNum = nrow(allSNPs)


#Output
write(paste("Total number of SNVs: ", totalNum, sep=""), file=paste(name, ".log", sep=""), append=T)
write(paste("After limiting by number of non-zeros: " , filteredNum, sep=""), file=paste(name, ".log", sep=""), append=T)


#Get size and labels of classes
if(is.null(args$sizeClass1)){
	
	stop("Did not provide size of first class.")
}else{

	num_class_1 = args$sizeClass1
	num_class_2 = ncol(allSNPs)-num_class_1
	class1Labels = colnames(allSNPs)[1:num_class_1]
	class2Labels = na.omit(colnames(allSNPs)[num_class_1+1:ncol(allSNPs)])
	write(paste("# of samples in ", args$nameClass1, ": ", num_class_1, sep=""), file=paste(name, ".log", sep=""), append=T)
	write(paste("Labels: ", paste(class1Labels, collapse=",", sep=""), sep=""), file=paste(name, ".log", sep=""), append=T)
	
	write(paste("# of samples in ", args$nameClass2, ": ", num_class_2, sep=""), file=paste(name, ".log", sep=""), append=T)
	write(paste("Labels: ", paste(class2Labels, collapse=",", sep=""), sep=""), file=paste(name, ".log", sep=""), append=T)
	classes = c(rep(1,num_class_1), rep(2, num_class_2))
	}
	



##Initializes the threads used for nested cross-validations.
cl<-makeCluster(args$threads)
registerDoParallel(cl)
write(paste("Number of threads:", args$threads, sep=" "), file=paste(name, ".log", sep=""), append=T)

#If provided, output to log which SNVs are being filtered
if(!(is.null(args$include))){
	write(paste("Including SNVs with substring:", args$include, sep=" "), file=paste(name, ".log", sep=""), append=T)
}
if(!(is.null(args$exclude))){
	write(paste("Excluding SNVs with substring: ", args$exclude, sep=""), file=paste(name, ".log", sep=""), append=T)
}


	
###If include or exclude is flagged, keep and remove those SNVs, respectively
if (    !(is.null(args$include))  |     !(is.null(args$exclude))   ){
        
	##If include is flagged...
	if(!(is.null(args$include))){
		#Splits the string by "," to identify all SNVs that are included
		types = strsplit(args$include, ",")[[1]]
		typed = c()
    for (i in 1:length(types)){
    	  #Finds the SNVs that have that substring
    		typed = union(rownames(allSNPs)[which(grepl(types[i],names[,2]))], typed)
    }

	}
	
	
	##If exlude is flagged..
	if (!(is.null(args$exclude))){
		#Splits the string by "," to identify all SNVs that are excluded
    excludes = strsplit(args$exclude, ",")[[1]]
		for (i in 1:length(excludes)){
          xor = rownames(allSNPs)[which(grepl(excludes[i],names[,2]))]
 			if(!(is.null(types[1])) ){
		  		typed = typed[which(!(typed %in% xor))]
 			}
			else{
					typed = allSNPs[which(!(rownames(allSNPs) %in% xor)),]
			}
		}
	
	}

}else{
	#If neither, set the found SNVs to be the total amount of SNVs
  typed = rownames(allSNPs)
}
allSNPs = allSNPs[which(rownames(allSNPs) %in% typed),]
write(paste("Number of SNVs in model (produced by filtering by type): ", length(typed), sep=""), file=paste(name, ".log", sep=""), append=T)
	

#Determine the number of CV iterations, create a matrix that contains the indicies of samples to be removed. One option is stratification
#where an equal number of samples are taken from each group. Another option is the amount of samples to take out for each iteration.

#If it is stratified...
if(!(args$unstratified)){

	#The total number of iterations is determined by this product
	maxIterations = choose(num_class_1, args$numOut) * choose(num_class_2, args$numOut)

	
	#If numCV is not specified, the number of iterations will be set to the maximum of maxIterations
	if(is.null(args$numCV)){
		numIter = maxIterations
		
	}
	else{
		#If numCV set, it set numIter to the minimum of the max and the specified value
		numIter = min(c(maxIterations, args$numCV))
	}
	
	#Produces a blank matrix fitting the dimensions of the study design and number of iterations
	skips = matrix(,nrow=numIter, ncol=args$numOut*2)
	if(numIter == maxIterations){
	
		index = 1
		#Gets every combination of the class labels for each class
		class1_combos = combn(seq(1,num_class_1,1), args$numOut)
		class2_combos = combn(seq(num_class_1+1,ncol(allSNPs),1), args$numOut)
		
		#For every combination of combinations...
		for (skip1 in 1:ncol(class1_combos)){
			for (skip2 in 1:ncol(class2_combos)){
				#Set the skip index matrix to that entry in the combo list
				skips[index,1:args$numOut] = class1_combos[1:args$numOut,skip1]
				skips[index,(args$numOut+1):(2*args$numOut)] = class2_combos[1:args$numOut,skip2]
				index = index + 1
			}
		}
	}else{
		#If numIterations is less than maximum iterations, we need to take a random subset of combinations.
		combos = c()
		
		for( i in 1:numIter){
			
			#Randomly choose combinations of class labels
			class1_index = sample(seq(1,num_class_1,1), args$numOut)
			class2_index = sample(seq(num_class_1+1, ncol(allSNPs),1),args$numOut)
			
			#Combine them togeter to get potentional set of indices
			combin = c(class1_index[1:args$numOut], class2_index[1:args$numOut])
			
			#This form loop goes through all previous combinations to ensure unique combinations
			found = F
			for(e in 1:length(combos)){
				numShared = length(union(combos[e], combin))
				if(numShared == 0){
					found=T
				}
			}
			#If the combo is unique, add it the set of combos
			if(!found){
				combos = cbind(combos, combin)
			}
		}
		#Convert the combos into a skip index matrix
	
		skips = t(data.frame(combos))
	
		}
	
}else{	#Samples should be unstratified

	maxIterations = choose(ncol(allSNPs), args$numOut) 
	
	#Set numIter to the minimum of the max iterations and user supplied amount
	if(is.null(args$numCV)){
		numIter = maxIterations
	}else{
		numIter = min(c(maxIterations, args$numCV))
	}
	#If the number of iterations is equal to the max amount of iterations..
	if(numIter == maxIterations){
		#Sets the skips index matrix to be all combinations of leave X out from all samples
		skips = t(combn(seq(1,ncol(allSNPs),1), args$numOut))
		}else{
		
		#Creates a list of unique combos chosen at random
		combos = c()
		while( length(combos) < numIter){
			combin = sample(seq(1,ncol(allSNPs),1), args$numOut)
			found = F
			for(e in 1:length(combos)){
				numShared = length(union(combos[e], combin))
				if(numShared == 0){
					found=T
				}
			}
			if(!found){
				combos = cbind(combos, combin)
			}
		}	
		#Sets the combos to be the skip index matrix
		skips = t(combos)
    }
}
write(paste("Maximum number of CV iterations: ", maxIterations, sep=""), file=paste(name, ".log", sep=""), append=T)
write(paste("Running number of CV iterations: ", numIter, sep=""), file=paste(name, ".log", sep=""), append=T)
if(numIter > 500){
	warning("Number of iterations is greater than 500. Runtime may be great.")
	write("Number of iterations is greater than 500. Runtime may be great.", file=paste(name, ".log", sep=""), append=T)	
}	
if(!args$unstrat){
	write(paste("Each CV removes ", args$numOut, " sample(s) from each class.", sep=""),file=paste(name, ".log", sep=""), append=T)
}else{
	write(paste("Each CV removes ", args$numOut, " sample(s) all samples.", sep=""),file=paste(name, ".log", sep=""), append=T)
}	


if(args$findOptimalK){
	
	write("Testing values of K.", file=paste(name, ".log", sep=""), append=T)
	#If a list of values of K are not provided
	if(is.null(args$range)){
		#Sees that the range provided makes sense
		
		testAmounts = tryCatch(seq(args$minK, args$maxK, args$everyKth), error=function(e){NULL})
		if(is.null(testAmounts)){
			
			stop("Range list not provided or range of K provide is erroneous.")	
		}
	}	else{
		subTestAmounts = tryCatch(seq(args$minK, args$maxK, args$everyKth), error=function(e){NULL})
		if(!is.null(subTestAmounts)){
			testAmounts = union(unlist(strsplit(args$range, ",")), subTestAmounts)
		}else{
			testAmounts = unlist(strsplit(args$range, ","))
		}
	}
	
	
	write(paste("Testing these values of K:", paste(testAmounts, sep=",",collapse=','), sep=" "), file=paste(name, ".log", sep=""), append=T)
	
	
	
	#bestK = c()
	bestK <- foreach (skip_index=1:nrow(skips), .packages=c('mixOmics', 'pROC')) %dopar% {
		#for(skip_index in 1:nrow(skips)){
		testIndex = skips[skip_index,]
		write(paste("Running Find Optimal K CV:",skip_index,sep=""), file=paste(name, ".log", sep=""), append=T)
		
		test.data = na.omit(allSNPs[,testIndex])
		if(length(testIndex) ==1){
			training.data = allSNPs[which(rownames(allSNPs) %in% attributes(test.data)$names),-testIndex]
			
		}else{
			training.data = allSNPs[which(rownames(allSNPs) %in% rownames(test.data)),-testIndex]
		}
		test.classes = classes[testIndex]
		training.data.classes = classes[-testIndex]
		sub.aucs = c()
		for (numFeatures in testAmounts){
			correct = 0
			missed = 0
			sub.class1 = c()
			resp = c()
			for (sub_skip in 1:ncol(training.data)){
				sub.test.class = training.data.classes[sub_skip]
				sub.training.classes = training.data.classes[-sub_skip]
				sub.test.data = na.omit(training.data[,sub_skip])
				keep = attributes(sub.test.data)$names
				sub.training.data = training.data[,-sub_skip][which(rownames(training.data) %in% keep),]
				#splsda.model <- tryCatch(splsda(t(sub.training.data), factor(t(sub.training.classes)), ncomp=1, keepX=numFeatures),  error = function(e) {NULL; print("fuck")})
				splsda.model <- splsda(t(sub.training.data), factor(t(sub.training.classes)), ncomp=1, keepX=as.numeric(numFeatures))
				if (is.null(splsda.model)){
					next
				}
				
				sub.test.data = sub.test.data[which(attributes(sub.test.data)$names %in% colnames(splsda.model$X))]
				prediction <- tryCatch(predict(splsda.model, t(sub.test.data)),  error = function(e) {NULL})
				if (is.null(prediction)){
					next
				}
				pred.classes <- prediction$class$centroids.dist;
				pred.class.one <- pred.classes[1];
				
				if (pred.class.one == sub.test.class){
					correct <- correct + 1;
				}else{
					missed <- missed + 1;
				}
				
				#Computing roc values
				unrounded.prediction.1 <- prediction$variates[1]
				class.1.centroid <- prediction$centroid[1]
				unrounded.prediction.1 <- unrounded.prediction.1 - class.1.centroid
				sub.class1 = c(sub.class1, unrounded.prediction.1)
				resp = c(resp, sub.test.class)
			}
			auc_val = as.numeric(ci.auc(roc(resp, sub.class1)))[2] + runif(1,-.0005,.0005)
			sub.aucs = c(sub.aucs , auc_val)       
		}
		
		subK = testAmounts[which(sub.aucs == max(sub.aucs))]
		subK
		#bestK = c(bestK, subK)
	}
	
	
	bestK = unlist(bestK)
	
	d = density(as.numeric(bestK))
	optK = round(d$x[which.max(d$y)])
	write(t(c(bestK)), ncol=1, file=paste(name, "_iteration_optimals.txt", sep=""), sep="\t")
  if(args$produceFigures){
  	pdf(paste(name,"_density.pdf",sep=""), 7,3)
  	dataset <- data.frame(X = bestK)
  	plot = ggplot(dataset, aes(x = X)) + geom_histogram(aes(y = ..density..), colour='black', fill='blue',binwidth=35) + geom_density(color='red', size=1)
  	plot = plot + geom_vline(xintercept=c(optM), linetype="dashed", size=1)
  	plot = plot + xlab("Value of K Tested")
  	plot = plot + ylab("Density")
  	plot = plot + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(limits = c(0,max(unlist(input))),expand = c(0,0))
  	plot = plot + theme_bw()
  	plot = plot + theme(axis.text=element_text(size=8), axis.title=element_text(size=10,face="bold"))
		plot
  	dev.off()	
  	}
	
	
	
	}else{
		if(is.null(args$optimalK)){
			stop("Did not specify opt. K or run --findOptK")
		}else{
			
			optK = as.numeric(args$optimalK)
		}
}
	

write(paste("Optimal K: ", optK, sep=""), file=paste(name, ".log", sep=""), append=T)


####Final CV using the best value of K
if(args$evaluatePerformance){
	correct <- 0
	missed <- 0
	OneWrong <- 0
	TwoWrong <- 0
	class1 <- c()
	class2 <- c()
	write("Evaluating Optimal K.", file=paste(name, ".log", sep=""), append=T)
	#performStats <- foreach (skip_index=1:nrow(skips), .packages=c('mixOmics', 'pROC')) %dopar% {
	for (skip_index in 1:nrow(skips)){   
		write(paste("Running Main CV:",skip_index,sep=""), file=paste(name, ".log", sep=""), append=T)
		testIndex = skips[skip_index,]
		test.data = na.omit(allSNPs[,testIndex])
		keep = rownames(test.data)
		if(length(testIndex) ==1){
			keep = attributes(test.data)$names
			
		}else{
			keep = rownames(test.data)
		}
		
		training.data = allSNPs[which(rownames(allSNPs) %in% keep),-testIndex]
		training.classes = classes[-testIndex]
		test.classes = classes[testIndex]
		splsda.model <- splsda(t(training.data), factor(t(training.classes)), ncomp=1, keepX=optK)
		
		training.data = training.data[which(rownames(training.data) %in% colnames(splsda.model$X)),]
		if(length(testIndex) ==1){
			test.data = test.data[which(attributes(test.data)$names %in% colnames(splsda.model$X))]
			
		}else{
			test.data = test.data[which(rownames(test.data) %in% colnames(splsda.model$X)),]
		}
		
		prediction <- predict(splsda.model, t(test.data))
		pred.classes <- prediction$class$centroids.dist;
		
		class.1.centroid <- prediction$centroid[1]
		class.2.centroid <- prediction$centroid[2]
		
		inv <- (class.1.centroid > class.2.centroid)
		
		if(inv){
			temp <- class.1.centroid
			class.1.centroid <- class.2.centroid
			class.2.centroid <- temp
		}
		
		class.2.centroid <- class.2.centroid - class.1.centroid
		
		unrounded.predictions.1 = c()
		unrounded.predictions.2 = c()
		
		for( i in 1:length(testIndex)){
			sampleClass = classes[testIndex[i]]
			if (pred.classes[i] == sampleClass){
				correct <- correct + 1;
				
				if(as.numeric(sampleClass) == 1){
					unrounded.predictions.1= c(unrounded.predictions.1, prediction$variates[i] - class.1.centroid / class.2.centroid)
				}else{
					unrounded.predictions.2= c(unrounded.predictions.2, prediction$variates[i] - class.1.centroid / class.2.centroid)
				}
				
			}else{
				missed <- missed + 1;
				if(as.numeric(sampleClass) == 1){
					
					OneWrong <- OneWrong + 1;
					unrounded.predictions.1= c(unrounded.predictions.1, prediction$variates[i] - class.1.centroid / class.2.centroid)
				}else{
					
					TwoWrong <- TwoWrong + 1	
					unrounded.predictions.2= c(unrounded.predictions.2, prediction$variates[i] - class.1.centroid / class.2.centroid)
				}
			}
		}
		
		if(inv){
			unrounded.predictions.1 = 1 - unrounded.predictions.1
			unrounded.predictions.2 = 1 - unrounded.predictions.2
		}
		
		class1_nim = "class1"
		class1_nim
		unrounded.predictions.1
		class2_nim = "class2"
		class2_nim
		unrounded.predictions.2
		class1 = c(class1, unrounded.predictions.1)
		class2 = c(class2, unrounded.predictions.2)
	}	
	
	#print(performStats)
	print(class1)
	print(class2)
	auc_ci = as.numeric(ci.auc(roc(controls=class1, cases=class2)))
	write(paste("AUC CI =  ", auc_ci, sep=""), file=paste(name, ".log", sep=""), append=T)
	
	accuracy =  correct / (correct + missed)
	print(accuracy)	
	write(paste("Accuracy: ", accuracy, sep=""), file=paste(name, ".log", sep=""), append=T)	
	
	stats = cbind.data.frame(c("Optimal K", "Predictive Accuracy", "Lower 95% CI","AUC","Upper 95% CI",paste(args$nameClass1, " Sensitivity", sep=""), paste(args$nameClass2, " Sensitivity", sep="")), c(optK, accuracy, auc_ci, (nrow(skips)-OneWrong)/nrow(skips), (nrow(skips)-TwoWrong)/nrow(skips))) 
	
	write(t(as.data.frame(stats)), ncol=2,file=paste(name, "_performance_stats.txt", sep=""), sep="\t")
	##Determine most relevant biomarkers with optimal M
	allData = allSNPs[which(rownames(allSNPs) %in% typed ),]
	total_classes = classes
	splsda.model <- splsda(t(allData), factor(t(total_classes)), ncomp=1, keepX=optk)
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
	
	topSNVs = t(snvs[order(as.numeric(snvs[,2]), decreasing = T),])
	write(topSNVs, ncol=3, file=paste(name, "_top_SNVs.txt", sep=""), sep="\t")
	
	
	if(args$produceFigures){
	  snv_names = c()
	  afs = c()
	  classOutcome = c()
	  for(e in 1:15){
	  	
	  	
	  	for(s in 1:ncol(allSNPs)){
	  	
	  		if(classes[s] == 1){
	  			classOutcome = c(classOutcome, nameClass1)	
	  			
	  		}else{
	  			classOutcome = c(classOutcome, nameClass2)
	  		}
	  		snv_name = topSNVs[e,1]
	  		snv_names = c(snv_names, snv_name) 
	  		afs = c(afs, allSNPs[snv_name,s])
	  	}
	  	
	  }
		
		
		df = data.frame(snp = factor(snv_names, levels=unique(snv_names)), af = afs, Outcome=classOutcome)
		
		pdf(paste(name,"_af_boxplot.pdf",sep=""), 7,3)
		plot4 = ggplot(df, aes(x=Outcome,y=af)) + geom_boxplot( frame.plot=FALSE,axes=FALSE,outlier.size=0) + facet_grid(. ~ snp) + facet_wrap( ~ snp, ncol=8)
		plot4 = plot4 +  geom_point(data=subset(df,af>0),size=1.25, aes(colour=Outcome),position = position_jitter(width = 0.1)) + scale_color_brewer(palette="Set1")
		plot4 = plot4 + theme_bw(base_size=10.5) 
		plot4 = plot4 + scale_y_continuous(breaks=c(0,.5,1))
		plot4 = plot4 + ylab("Allele Fraction")
		plot4 = plot4 + xlab(NULL) + ggtitle(paste("Top ", name, " Nonsynonymous Exonic SNVs",sep=""))
		plot4 = plot4 + theme(axis.title.x=element_blank(),
													axis.text.x=element_blank(),
													axis.ticks.x=element_blank(), legend.position="bottom")
		plot4
		dev.off()
	
		pdf(paste(name,"_heatmap.pdf", sep=""), 7,5)
		data = allSNPs[which(rownames(allSNPs) %in% topSNVs[1:15,]),]
		SNVs = rownames(allSNPs)
		samps = colnames(allSNPs)
		colnames(data) = samps
		rownames(data) = t(SNVs)
		data = matrix(as.numeric(unlist(data)), nrow=nrow(data))
		colnames(data) = samps
		rownames(data) = t(SNVs)
		heatmap.2(data, key.title="NA", key.xlab="Allele Fraction",density.info='none',na.color='black', col = colorpanel(100,"white","blue"),margins=c(4.2,10.4),trace='none')
		dev.off()
		}
}

if(args$permutationTests){

	if(!(args$findOptimalK) & is.null(args$optimalK)){
		stop("Need to provide optimal K or run findOptimalK")	
	}
	if(!(args$evaluatePerformance) & is.null(args$AreaUnderCurve)){
		stop("Need to provide AUC or run evaluatePerformance")	
	}
	if(!is.null(args$AreaUnderCurve)){
	   compare_auc = args$AreaUnderCurve	
	}else{
		compare_auc = auc_ci[2]
	}

	numIter = as.numeric(args$permutationIterations)
	write(paste("Running ",numIter," permuation tests.",sep=""),file=paste(name, ".log", sep=""), append=T)	

	AUCs = c()
	run = 0
	ls <- foreach (iter=1:numIter, .packages=c('mixOmics','pROC')) %dopar% {
		#for(iter in 1:numIter){
		write(paste("Iteration: ",iter, sep=""), file=paste(name, ".log", sep=""), append=T)
		
		correct <- 0
		missed <- 0
		OneWrong <- 0
		TwoWrong <- 0
		class1 <- c()
		class2 <- c()
		allSNPs = allSNPs[, sample(c(seq(1,ncol(allSNPs),1)))]
		
		#performStats <- foreach (skip_index=1:nrow(skips), .packages=c('mixOmics', 'pROC')) %dopar% {
		for (skip_index in 1:nrow(skips)){   
			testIndex = skips[skip_index,]
			test.data = na.omit(allSNPs[,testIndex])
			keep = rownames(test.data)
			if(length(testIndex) ==1){
				keep = attributes(test.data)$names
				
			}else{
				keep = rownames(test.data)
			}
			
			training.data = allSNPs[which(rownames(allSNPs) %in% keep),-testIndex]
			training.classes = classes[-testIndex]
			test.classes = classes[testIndex]
			splsda.model <- splsda(t(training.data), factor(t(training.classes)), ncomp=1, keepX=optK)
			
			training.data = training.data[which(rownames(training.data) %in% colnames(splsda.model$X)),]
			if(length(testIndex) ==1){
				test.data = test.data[which(attributes(test.data)$names %in% colnames(splsda.model$X))]
				
			}else{
				test.data = test.data[which(rownames(test.data) %in% colnames(splsda.model$X)),]
			}
			
			prediction <- predict(splsda.model, t(test.data))
			pred.classes <- prediction$class$centroids.dist;
			
			class.1.centroid <- prediction$centroid[1]
			class.2.centroid <- prediction$centroid[2]
			
			inv <- (class.1.centroid > class.2.centroid)
			
			if(inv){
				temp <- class.1.centroid
				class.1.centroid <- class.2.centroid
				class.2.centroid <- temp
			}
			
			class.2.centroid <- class.2.centroid - class.1.centroid
			
			unrounded.predictions.1 = c()
			unrounded.predictions.2 = c()
			
			for( i in 1:length(testIndex)){
				sampleClass = classes[testIndex[i]]
				if (pred.classes[i] == sampleClass){
					correct <- correct + 1;
					
					if(as.numeric(sampleClass) == 1){
						unrounded.predictions.1= c(unrounded.predictions.1, prediction$variates[i] - class.1.centroid / class.2.centroid)
					}else{
						unrounded.predictions.2= c(unrounded.predictions.2, prediction$variates[i] - class.1.centroid / class.2.centroid)
					}
					
				}else{
					missed <- missed + 1;
					if(as.numeric(sampleClass) == 1){
						
						OneWrong <- OneWrong + 1;
						unrounded.predictions.1= c(unrounded.predictions.1, prediction$variates[i] - class.1.centroid / class.2.centroid)
					}else{
						
						TwoWrong <- TwoWrong + 1	
						unrounded.predictions.2= c(unrounded.predictions.2, prediction$variates[i] - class.1.centroid / class.2.centroid)
					}
				}
			}
			
			if(inv){
				unrounded.predictions.1 = 1 - unrounded.predictions.1
				unrounded.predictions.2 = 1 - unrounded.predictions.2
			}
			
			class1_nim = "class1"
			class1_nim
			unrounded.predictions.1
			class2_nim = "class2"
			class2_nim
			unrounded.predictions.2
			class1 = c(class1, unrounded.predictions.1)
			class2 = c(class2, unrounded.predictions.2)
		}	
		auc_ci = as.numeric(ci.auc(roc(controls=class1, cases=class2)))
		auc_ci[2]
	}

aucs = unlist(ls)
write(unlist(ls), ncol=1, file=paste(name, "_permutation.txt", sep=""), sep="\t")
p_val = length(which(aucs >= compare_auc))/length(aucs)
write(paste("Perm p-value: ",p_val, sep=""), file=paste(name, ".log", sep=""), append=T)
write(p_val, ncol=1, file=paste(name,"_p_val.txt", sep=""),sep="\t")

}
stopCluster(cl)		
	




