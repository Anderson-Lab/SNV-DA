suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))



help_text = "\nSNV-DA is used to create and evaluate sPLS-DA models to identify single nucleotide variations (SNVs) that accurately classify phenotype. The steps of the pipeline include:\n1) Finding optimal value of number of features to be selected in the model (K)\n2) Evaluating the model using cross-validations and optimal value of K\n3) Find and rank the K features that are correlated with predictive accuracy\n4) Permutation tests to determine if model is discriminate towards the true grouping of labels.\n\n
The script has several parameters that allows the user to design their own analysis:\n1) The number of non-zero values by which to filter the matrix\n2) The number of NA values allowed for each SNV\n3) The type of SNV to include in the model\n4) The range and/or values of K tested.\n4) Cross-validation design (number of samples to take out from each group, the stratification of test samples [force test samples to be pulled equally from each group], and the number of cross-validations).\n\nHappy hunting!"  
option_list <- list(
make_option(c("-O", "--findOptimalK"), action="store_true", default=F,
								 help="If flagged, SNV-DA will run nest cross-validations to determine optimal K."),

make_option(c("-J", "--evaluatePerformance"), action="store_true", default=F,
								 help="If flagged, SNV-DA will evaluate performance of model the value of K found from running --findOptimalK or user supplied --optimalK."),

make_option(c("-P", "--permutationTests"), action="store_true", default=F,
								 help="If flagged, SNV-DA will run permutation tests by randomly permuting sample classes then running cross-validations. Tests are run --permIter number of times. Performances are then compared to the value of AUC found by --evaluatePerformance of user supplied --AreaUnderCurve."),

make_option(c("-x", "--predictTestSet"), action="store_true", default=F,
						help="If flagged, SNV-DA will train a model with optimal K on all samples in the main matrix and then predict the labels of the second matrix."),


make_option(c("-C", "--getSNVs"), action="store_true", default=F,
						help="If flagged, SNV-DA will ONLY run sPLS-DA on all samples to identify K top predictive features"),


make_option(c("-W", "--NfoldCV"), type="integer", default=NULL,
                      help="If specified, N-fold cross-validations will be performed by partitioning the input matrix into NfoldCV setsi during the optimization of K and performance evaluation."),

make_option(c("-I", "--iterNfoldCV"), type="integer", default=1,
                      help="If NfoldCV is specified, NfoldCV will be run iterNFoldCV times during performance evaluation  and finding of optimal K (default=1). "),
make_option(c("-V", "--internalNfoldCV"), type="integer", default=NULL,
						help="If specified, N-fold cross-validations will be performed during the interval CV while finding optimal K by partitioning the input matrix into iternalNfoldCV setsi during the optimization of K and performance evaluation. If left unspecified, leave-one-out crossvalidations will be used for the interval cross-validation."),

make_option(c("-j", "--internalIterNfoldCV"), type="integer", default=1,
						help="If internalNfoldCV is specified, internalNfoldCV will be run internalIterNFoldCV times during performance evaluation and finding of optimal K (default=1)."),

make_option(c("-M","--matrix"),type="character", default=NULL,
								 help = "SNVM input (csv)"),
make_option(c("-y","--testMatrix"),type="character", default=NULL,
 						help = "matrix to predict labels (csv)"),

make_option(c("-D","--sizeClass1"),type="integer", default=NULL,
								 help = "The number of samples of class 1. (Used to find the number of samples of both classes)."),

make_option(c("-L", "--minK"), type="integer", default=NULL,
								 help="The minimum value of K tested."),

make_option(c("-B", "--maxK"), type="integer",default=NULL,
								 help="The maximum value of K tested."),

make_option(c("-E", "--everyKth"), type="integer",default=NULL,
								 help="Specifies every Nth K tested.\nOr use:"),

make_option(c("-R", "--range"), type="character",default=NULL,
								 help="A comma separated list of K's to be added to the Ks tested, e.g. 5,10,15,20,25,50,75,100,250,300."),

make_option(c("-Z", "--numNonZeros"), type="integer", default=3,
								 help="Limit SNVs to those that have at least X number of samples with allele fractions greater than zero (default = 3). Decreases runtime. Suggested value is a 1/3 the number of samples in the smaller group."),

make_option(c("-N", "--numNAs"), type="integer", default=NULL,
								 help="Limit SNVs to those that have at most X amount of samples with 'NA' values (default=0)."),

make_option(c("-n", "--NonNormImpute"), action="store_true", default=F,
						help="If flagged, NAs are imputed with non-standardized values."),


make_option(c("-G", "--include"), type="character", default=NULL,
								 help="SNVs with this substring are included, e.g. \"exonic\" for nonsynonymous and synonymous exonic SNVs"),

make_option(c("-X", "--exclude"), type="character", default=NULL,
								 help="SNVs with this substring are excluded, e.g. ncRNA to remove intronic_ncRNA"),

make_option(c("-A","--studyName"),type="character", default="SNV-DA_run",
						help = "One word name of the analysis used in output file names and figures."),

make_option(c("-U","--nameClass1"),type="character", default="Class_1",
						help = "Name of class 1 (default = Class_1)."),

make_option(c("-H","--nameClass2"),type="character", default="Class 2",
						help = "Name of class 2 (default = Class_2)."),

make_option(c("-T", "--threads"), type="integer", default=1,
						help="Number of threads for cross-validations (default = 1)."),

make_option(c("-K", "--optimalK"), type="integer", default=500,
						help="Used when --evaluatePerformance or --permuationTests is flagged. Specifies the value of K used for cross-validation (default = 500)."),

make_option(c("-q", "--permIter"), type="integer", default=NULL,
						help="Used when --permuationTests is flagged. Specifies the number of permutation tests to be run"),

make_option(c("-a", "--AreaUnderCurve"), type="double", default=NULL,
						help="Performance of true model. Compared to permuted models when --permutationTests are run. If unflagged, value is determined by --evaluatePerformance"),

make_option(c("-F", "--produceFigures"), action="store_true", default=F,
						help="If flagged, allele fraction box plots, heatmaps, and kernal density figures are produced.")

)

args <- parse_args(OptionParser(description=help_text, option=option_list))

get_predictions <- function(skips,training,tmp.classes, K_hat){
	
	correct <- 0
	missed <- 0
	OneWrong <- 0
	TwoWrong <- 0
	numOne = 0
	numTwo = 0
	class_hat <- c()
	class  = c()

	for (skip_index in 1:nrow(skips)){   

		instance = skips[skip_index,]
	  instance = instance[!is.na(instance)]
		testIndex = as.numeric(instance)
		test.data = data.frame(training[,testIndex])


		rownames(test.data) = rownames(training)

		training.data = training[,-testIndex]
		training.classes = tmp.classes[-testIndex]
		test.classes = classes[testIndex]
		
		num_class_1 = length(training.classes[which(training.classes == 1)])
		num_class_2 = length(training.classes) - num_class_1

		training.data.1 = training.data[,1:num_class_1]
		training.data.2 = training.data[,(num_class_1+1):ncol(training.data)]

		###Removes any features where either group has all NAs (non-informative features)
		training.data = training.data[which(rowSums(is.na(training.data.1)) < num_class_1 & rowSums(is.na(training.data.2)) < num_class_2),]
		###Remove any features where all AFs are the same (non-informatic features)
		training.data = training.data[which(apply(training.data,1, function(x) length(na.omit(unique(x)))) > 1),]
		
		test.data = test.data[which(rowSums(is.na(test.data)) != ncol(test.data)),]
		
		training.data = training.data[which(rownames(training.data) %in% rownames(test.data)),]

		splsda.model <- tryCatch(splsda(t(training.data), factor(t(training.classes)), ncomp=1, keepX=K_hat),  error = function(e) {NULL})
 		if (is.null(splsda.model)){
 			next
		}

		#Only keep test.data of features that are still included in model
		test.data = as.data.frame(test.data[which(rownames(test.data) %in% colnames(splsda.model$X)),])
  	training.data = as.data.frame(training.data[which(rownames(training.data) %in% colnames(splsda.model$X)),])
		
  	#Use non-normalized values for imputation
  	if(args$NonNormImpute){
		##Get the training data that correspond to features that NA values in the test set
		class1_na_rows = training.data[,1:num_class_1]
		class2_na_rows = training.data[,(num_class_1+1):ncol(training.data)]

  	}else{
  		
  		class1_na_rows = t(splsda.model$X[1:num_class_1,])
  		class1_na_rows[which(apply(class1_na_rows, 1, function(x) all(is.na(x)))),] = rep(0, ncol(class1_na_rows))
  		class2_na_rows = t(splsda.model$X[(num_class_1+1):ncol(training.data),])
  		class2_na_rows[which(apply(class2_na_rows, 1, function(x) all(is.na(x)))),] = rep(0, ncol(class2_na_rows))
  		
  	}
  	##Get mean of means for imputations
  	imputes = rowMeans(cbind(rowMeans(class1_na_rows, na.rm=T), rowMeans(class2_na_rows, na.rm=T)), na.rm=T)
		
  	test.data.impute = test.data

		##For each sample in the test set, fill the NA values with imputed values from the training set (mean of means)
		for(samp in 1:ncol(test.data)){
			nas = which(is.na(test.data[,samp]))
			test.data.impute[nas, samp] = as.numeric(imputes[nas])
		}


		
 		prediction <- tryCatch(predict(splsda.model, t(test.data.impute)),  error = function(e) {NULL})
 		if (is.null(prediction)){
 			next
 		}
		
		pred.classes <- prediction$class$centroids.dist
		class.1.centroid <- prediction$centroid[1]
		class.2.centroid <- prediction$centroid[2]
		
		unrounded.predictions = c()
		
		for( i in 1:length(testIndex)){
			sampleClass = test.classes[i]
			if (pred.classes[i] == sampleClass){
				correct <- correct + 1;
			}else{
				missed <- missed + 1;
				if(as.numeric(sampleClass) == 1){
					OneWrong <- OneWrong + 1;
				}else{
					TwoWrong <- TwoWrong +1
				}
			}
			if (sampleClass == 1){
				numOne = numOne +1	
			}else{
				numTwo = numTwo +1	
			}
			
			unrounded.prediction <- (prediction$variates[i] - class.1.centroid) / (class.2.centroid - class.1.centroid)
			mx_label = max(training.classes)
			mn_label = min(training.classes)
			unrounded.prediction = unrounded.prediction * (mx_label-mn_label) + mn_label
			unrounded.predictions = c(unrounded.predictions, unrounded.prediction)
		}
		
		
		class_hat = c(class_hat, unrounded.predictions)			
		class = c(class,  test.classes)
	}	
	

	AUC = as.numeric(ci.auc(roc(class, class_hat)))
	ACC = correct / (correct + missed)
	Sens.1 = (numOne-OneWrong)/numOne
	Sens.2 = (numTwo-TwoWrong)/numTwo
	
	list(AUC = AUC, ACC = ACC, Sens.1 = Sens.1, Sens.2 = Sens.2)

}


#Argument handling
name = args$studyName

##Check to see if SNVM was input
if(is.null(args$matrix)){
	
stop("SNV matrix not provided")	
}
write(paste("Running analysis: ", name, sep=" "), file=paste(name, ".log", sep=""), append=T)
write(paste("Input SNVM:", args$matrix, sep=" "), file=paste(name, ".log", sep=""), append=T)

##Output to the log the required number of nonzero values allows
write(paste("Removing SNVs that do not have at least ", args$numNonZeros," samples with allele fraction values >0", sep=""), file=paste(name, ".log", sep=""), append=T)
if(!is.null(args$numNAs)){
	write(paste("Removing SNVs that have more than ", args$numNAs," samples with NA values.", sep=""), file=paste(name, ".log", sep=""), append=T)
}


##Read in SNVM and set the names of the loci. 
allSNPs <- read.csv(args$matrix, header=T, stringsAsFactors=F,sep=",",  na.strings = "Na")
names <- allSNPs[,c(1,2)]
rownames(names) = names[,1]
allSNPs <- data.matrix(allSNPs[,c(-1,-2)])
rownames(allSNPs) = names[,1]
totalNum = nrow(allSNPs)

##Determine rows that have X non-zero feature values

if(is.null(args$numNAs)){
	numNAs = ncol(allSNPs)
}else{numNAs = args$numNAs}

if(args$numNonZeros > 1 | numNAs < ncol(allSNPs)){
	
	indices<-which(rowSums(is.na(allSNPs))<=numNAs & rowSums(allSNPs>0,TRUE)>=args$numNonZeros)
}else{indices = rownames(allSNPs)}


##Pull out those SNVs
names = names[indices,]
allSNPs = allSNPs[indices,]
filteredNum = nrow(allSNPs)

#Output
write(paste("Total number of SNVs: ", totalNum, sep=""), file=paste(name, ".log", sep=""), append=T)
write(paste("After limiting by number of non-zeros and NA values: " , filteredNum, sep=""), file=paste(name, ".log", sep=""), append=T)

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
registerDoParallel(cores=args$threads)
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
	if(!(is.null(c))){
		write(paste("Only including SNVs whose type has substring: ", args$include, sep=""), file=paste(name, ".log", sep=""), append=T)
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
		write(paste("Excluding SNVs whose type has substring: ", args$exclude, sep=""), file=paste(name, ".log", sep=""), append=T)
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
	

#Determine the number of CV iterations, create a matrix that contains the indicies of samples to be removed.
if(args$predictTestSet){
	
	training.data = allSNPs[which(rownames(allSNPs) %in% typed ),]
	total_classes = classes
	print(paste("Predicting classes of test set with optimal K =", args$optimalK))
	num_class_1 = args$sizeClass1
	num_class_2 = ncol(allSNPs)-num_class_1
	training.data.1 = training.data[,1:num_class_1]
	training.data.2 = training.data[,(num_class_1+1):ncol(training.data)]
	
	training.data = training.data[which(rowSums(is.na(training.data.1)) < num_class_1 & rowSums(is.na(training.data.2)) < num_class_2),]
	###Remove any features where all AFs are the same (non-informatic features)
	training.data = training.data[which(apply(training.data,1, function(x) length(na.omit(unique(x)))) > 1),]
	
	testSNPs <- read.csv(args$testMatrix, header=T, stringsAsFactors=F,sep=",",  na.strings = "Na")
	names <- testSNPs[,c(1,2)]
	rownames(names) = names[,1]
	suppressWarnings(testSNPs <- data.matrix(testSNPs[,c(-1,-2)]))
	rownames(testSNPs) = names[,1]

	testSNPs = testSNPs[which(rowSums(is.na(testSNPs)) <= args$numNAs),]

	training.data = training.data[which(rownames(training.data) %in% rownames(testSNPs)),]
	
	
	splsda.model <- splsda(t(training.data), factor(t(total_classes)), ncomp=1, keepX=args$optimalK)
	SNVs = splsda.model$loadings$X
	all_snv_names = rownames(SNVs)	
	training.data = training.data[all_snv_names,]
	
	subtestSNPs = testSNPs[all_snv_names,]

	numTest1 = sum(grepl(gsub("-", ".", args$nameClass1), colnames(subtestSNPs)))
	numTest2 = length(colnames(subtestSNPs)) - numTest1
	test.classes = c(rep(1, numTest1), rep(2, numTest2))

	class1_na_rows = training.data[,1:num_class_1]
	class1_na_rows[which(apply(class1_na_rows, 1, function(x) all(is.na(x)))),] = rep(0, ncol(class1_na_rows))
	
	class2_na_rows = training.data[,(num_class_1+1):ncol(training.data)]
	class2_na_rows[which(apply(class2_na_rows, 1, function(x) all(is.na(x)))),] = rep(0, ncol(class2_na_rows))
	
	##Get mean of means for imputations
	imputes = rowMeans(cbind(rowMeans(class1_na_rows, na.rm=T), rowMeans(class2_na_rows, na.rm=T)))
 	
	##For each sample in the test set, fill the NA values with imputed values from the training set (mean of means)

	for(samp in 1:ncol(subtestSNPs)){
		nas = which(is.na(subtestSNPs[,samp]))
		subtestSNPs[nas, samp] = as.numeric(imputes[nas])
	}

	prediction = predict(splsda.model, t(subtestSNPs))
	pred.classes <- prediction$class$centroids.dist
	output = t(cbind(colnames(testSNPs), pred.classes))

	write.table(output, file=paste(args$studyName,"_",args$include,"_ext_pred_results.txt", sep=""), col.names=F, row.names=F,quote=F, append=T)
	class.1.centroid <- prediction$centroid[1]
	class.2.centroid <- prediction$centroid[2]
	
	unrounded.predictions = c()
	correct = 0
	missed = 0
	OneWrong = 0
	TwoWrong = 0
	numOne = 0
	numTwo = 0
	print(test.classes)
	print(pred.classes)
	for( i in 1:length(test.classes)){
		sampleClass = test.classes[i]
		if (pred.classes[i] == sampleClass){
			correct <- correct + 1;
		}else{
			missed <- missed + 1;
			if(as.numeric(sampleClass) == 1){
				OneWrong <- OneWrong + 1;
			}else{
				TwoWrong <- TwoWrong +1
			}
		}
		if (sampleClass == 1){
			numOne = numOne +1	
		}else{
			numTwo = numTwo +1	
		}
		
		
		unrounded.prediction <- (prediction$variates[i] - class.1.centroid) / (class.2.centroid - class.1.centroid)
		mx_label = max(total_classes)
		mn_label = min(total_classes)
		unrounded.prediction = unrounded.prediction * (mx_label-mn_label) + mn_label
		unrounded.predictions = c(unrounded.predictions, unrounded.prediction)	
	}
	print(unrounded.predictions)
	AUC = as.numeric(ci.auc(roc(test.classes, unrounded.predictions)))
	ACC = correct / (correct + missed)
	Sens.1 = (numOne-OneWrong)/numOne
	Sens.2 = (numTwo-TwoWrong)/numTwo
	write(paste("testset AUC: ", AUC), file=paste(args$studyName,"_",args$include,"_ext_pred_results.txt", sep=""), append=T)
	write(paste("testset predictive accuracy: ", ACC), file=paste(args$studyName,"_",args$include,"_ext_pred_results.txt", sep=""), append=T)
	write(paste("testset class 1 accuracy: ", Sens.1), file=paste(args$studyName,"_",args$include,"_ext_pred_results.txt", sep=""),  append=T)
	write(paste("testset class 2 accuracy: ", Sens.2), file=paste(args$studyName,"_",args$include,"_ext_pred_results.txt", sep=""),  append=T)
	
	print(paste("testset AUC: ", AUC))
	print(paste("testset predictive accuracy: ", ACC))
	print(paste("testset class 1 accuracy: ", Sens.1))
	print(paste("testset class 2 accuracy: ", Sens.2))
	
}
if(!args$predictTestSet){
if(is.null(args$NfoldCV)){

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
}else{

	#Setting up unstratified N-fold cross-validations and storing the indexes of a sample in the skips matrix
	indices = 1:ncol(allSNPs)
	numPer = floor(ncol(allSNPs)/args$NfoldCV)
	rem = ncol(allSNPs) %% args$NfoldCV
	numIter = args$iterNfoldCV
	skips = matrix(NA, nrow=args$NfoldCV*numIter, ncol=numPer+1)
	for(it in 0:numIter-1){
		rando = sample(1:ncol(allSNPs))
		for(nf in 1:args$NfoldCV){
			if(nf <= rem){
				
				skips[((it*args$NfoldCV)+nf),] = rando[((numPer+1)*nf-numPer):(nf*(numPer+1))]
				
			}else{
				skips[((it*args$NfoldCV)+nf),] = c(rando[(numPer*nf-numPer+1+rem):(nf*(numPer)+rem)], NA)
			}
		}
	}
	skips = data.frame(skips)
	
}
		
}

# write(paste("Maximum number of CV iterations: ", maxIterations, sep=""), file=paste(name, ".log", sep=""), append=T)


if(args$findOptimalK && !args$predictTestSet){

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
	

	bestK <- foreach (skip_index=1:nrow(skips), .packages=c('mixOmics', 'pROC')) %dopar% {
	#bestK = c()
	#for(skip_index in 1:nrow(skips)){
	 	
	  skips2 = skips[skip_index,]
	  skips2 = skips2[!is.na(skips2)]
		testIndex = skips2
		write(paste("Running Find Optimal K CV:",skip_index,sep=""), file=paste(name, ".log", sep=""), append=T)
		test.data = allSNPs[,testIndex]
		training.data = allSNPs[,-testIndex]
		test.classes = classes[testIndex]
		training.data.classes = classes[-testIndex]
		sub.aucs = c()

		for (numFeatures in testAmounts){

			write(paste("CV: ",skip_index, " Evaluating K: ", numFeatures,sep=""), file=paste(name, "_progress.log", sep=""), append=T)
			if(is.null(args$internalNfoldCV)){
			
				results = get_predictions(as.data.frame(c(1:ncol(training.data))), training.data, training.data.classes, numFeatures)
				sub.aucs = c(sub.aucs , results$AUC[2])       
			}else{
			
					#Setting up unstratified NfoldCV for internal CV
					indices = 1:ncol(training.data)
					numPer = floor(ncol(training.data)/args$internalNfoldCV)
				
					rem = ncol(training.data) %% args$internalNfoldCV
					internal.skips = matrix(NA, nrow=args$internalNfoldCV, ncol=numPer+1)
					it = 0
					rando = sample(1:ncol(training.data))
					for(nf in 1:args$internalNfoldCV){
						if(nf <= rem){
							
							internal.skips[((it*args$internalNfoldCV)+nf),] = rando[((numPer+1)*nf-numPer):(nf*(numPer+1))]
							
						}else{
							internal.skips[((it*args$internalNfoldCV)+nf),] = c(rando[(numPer*nf-numPer+1+rem):(nf*(numPer)+rem)], NA)
						}
					}
					internal.skips = data.frame(internal.skips)
						

				results = get_predictions(internal.skips, training.data, training.data.classes, numFeatures)
				sub.aucs = c(sub.aucs , results$AUC[2])  
				
			  }
			}
		
		subK = testAmounts[which(sub.aucs == max(sub.aucs))]
		write(subK, file=paste(name, "_optK_iteration_optimals.txt", sep=""), sep="\t", append=T)
 		#bestK = c(bestK, subK)
 		subK
	}
	
	bestK = c(as.numeric(unlist(bestK)))
	d = density(as.numeric(bestK))
	optK = round(d$x[which.max(d$y)])
	write(t(c(bestK)), ncol=1, file=paste(name, "_iteration_optimals.txt", sep=""), sep="\t")
  if(args$produceFigures){
  	pdf(paste(name,"_density.pdf",sep=""), 7,3)
  	dataset <- data.frame(X = as.numeric(bestK))
    if(is.null(args$everyK)){
    	binwidth=args$everyK
    	}else{
    		binwidth=5
    	}
  	
  	plot = ggplot(dataset, aes(x = X)) + geom_histogram(aes(y = ..density..), colour='black', fill='blue',binwidth=binwidth) + geom_density(color='red', size=1)
  	plot = plot + geom_vline(xintercept=c(optK), linetype="dashed", size=1)
  	plot = plot + xlab("Value of K Tested")
  	plot = plot + ylab("Density")
  	plot = plot + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(limits = c(0,max(as.numeric(bestK))),expand = c(0,0))
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
if(args$evaluatePerformance && !args$predictTestSet){
	if(args$iterNfoldCV == 1){
		
		results = get_predictions(skips, allSNPs, classes, optK)
	
		auc_ci = results$AUC
		acc = results$ACC
		class1.sens = results$Sens.1
		class2.sens = results$Sens.2
	}else{
	
		aucs = c()
		acc = c()
		class1.sens = c()
		class2.sens = c()
		eval <- foreach (iter=0:(args$iterNfoldCV-1), .packages=c('mixOmics', 'pROC')) %dopar% {
		#for(iter in 0:(args$iterNfoldCV-1)){
				write(paste("Evaluating model, iteration: ",(iter+1),sep=""), file=paste(name, ".log", sep=""), append=T)
				iter_skips = skips[(iter*args$NfoldCV + 1):(iter*args$NfoldCV + args$NfoldCV),]
				results = get_predictions(iter_skips, allSNPs, classes, optK)
				
				#aucs = c(aucs, results$AUC[2])
				write(results$AUC[2], file=paste(name, "_AUCs.log", sep=""), append=T)
				write(c(results$ACC, results$Sens.1, results$Sens.2), sep="\t",file=paste(name, "_preds.log", sep=""), append=T)
				#acc = c(acc, results$ACC)
				#class1.sens = c(class1.sens, results$Sens.1)
				#class2.sens = c(class2.sens, results$Sens.2)
				list(results$AUC[2],results$ACC, results$Sens.1, results$Sens.2)
		}
		

		eval = matrix(unlist(eval), ncol=args$iterNfoldCV)
		
		AUC = mean(eval[1,])
		std.err = sd(eval[1,])/sqrt(length(unlist(eval[1,])))
		
		auc_ci = c(AUC-1.96*std.err, AUC, AUC+1.96*std.err)
		acc = mean(eval[2,])
		class1.sens = mean(eval[3,])
		class2.sens = mean(eval[4,])

	}	

	write(paste("AUC CI =  ", auc_ci, sep=""), file=paste(name, ".log", sep=""), append=T)
	write(paste("Accuracy: ", acc, sep=""), file=paste(name, ".log", sep=""), append=T)	
	
	stats = cbind.data.frame(c("Optimal K", "Predictive Accuracy", "Lower 95% CI","AUC","Upper 95% CI",paste(args$nameClass1, " Sensitivity", sep=""), paste(args$nameClass2, " Sensitivity", sep="")), c(optK, acc, auc_ci, class1.sens, class2.sens)) 
	
	write(t(as.data.frame(stats)), ncol=2,file=paste(name, "_performance_stats.txt", sep=""), sep="\t")
	##Determine most relevant biomarkers with optimal M
	allData = allSNPs[which(rownames(allSNPs) %in% typed ),]
	total_classes = classes
	splsda.model <- splsda(t(allData), factor(t(total_classes)), ncomp=1, keepX=optK)
	SNVs = splsda.model$loadings$X
	all_snv_names = rownames(SNVs)
	snvs = matrix(, nrow=optK, ncol=3)
	index = 1
	for(i in 1:length(SNVs)){
		if (abs(SNVs[i]) > 0){
			snvs[index,1:2] = c(all_snv_names[i], SNVs[i])
			index = index + 1
		}
	}
	snvs[,2] = as.numeric(snvs[,2])
	topSNVs = snvs[order(abs(as.numeric(snvs[,2])), decreasing = T),]
	
	rownames(topSNVs) = topSNVs[,1]
	topSNVs[,3] = names[rownames(topSNVs),2]
	write(t(topSNVs), ncol=3, file=paste(name, "_top_SNVs.txt", sep=""), sep="\t")
	write(t(topSNVs), ncol=3, file=paste(name, "_top_SNVs.csv", sep=""), sep=",")
	
	
	if(args$produceFigures){
	  snv_names = c()
	  afs = c()
	  classOutcome = c()
	  numTop = min(optK, 15)
	  for(e in 1:numTop){
	  	for(s in 1:ncol(allSNPs)){
	  		
	  		snv_name = topSNVs[e,1]
	  		af = allSNPs[snv_name,s]
	  		if(classes[s] == 1){
	  			classOutcome = c(classOutcome, args$nameClass1)	
	  		}else{
	  			classOutcome = c(classOutcome, args$nameClass2)
	  		}
	  		snv_names = c(snv_names, snv_name) 
	  		afs = c(afs, af)
	  		
	  	}
	  }
		
		df = data.frame(snp = factor(snv_names, levels=unique(snv_names)), af = afs, Outcome=classOutcome)
		df = na.omit(df)
		pdf(paste(name,"_af_boxplot.pdf",sep=""), 7,3)
		plot4 = ggplot(df, aes(x=Outcome,y=af)) + geom_boxplot( frame.plot=FALSE,axes=FALSE,outlier.size=0) + facet_grid(. ~ snp) + facet_wrap( ~ snp, ncol=min(ceiling(numTop/2), 8))
		plot4 = plot4 +  geom_point(data=subset(df,af>0),size=1.25, aes(colour=Outcome),position = position_jitter(width = 0.1)) + scale_color_brewer(palette="Set1")
		plot4 = plot4 + theme_bw(base_size=5.5) 
		plot4 = plot4 + scale_y_continuous(breaks=c(0,.5,1))
		plot4 = plot4 + ylab("Allele Fraction")
		plot4 = plot4 + xlab(NULL) + ggtitle(paste("Top ", name, " SNVs",sep=""))
		plot4 = plot4 + theme(axis.title.x=element_blank(),
													axis.text.x=element_blank(),
													axis.ticks.x=element_blank(), legend.position="bottom")
		plot4
		print(plot4)
		dev.off()
	
		pdf(paste(name,"_heatmap.pdf", sep=""), 8,5)
		data = allSNPs[which(rownames(allSNPs) %in% topSNVs[1:numTop,]),]
		SNVs = rownames(data)
		samps = colnames(allSNPs)
		data = matrix(as.numeric(unlist(data)), nrow=nrow(data))
		colnames(data) = samps
		rownames(data) = t(SNVs)
		heatmap.2(data, key.title="NA", key.xlab="Allele Fraction",density.info='none',na.color='black', col = colorpanel(100,"white","blue"),margins=c(4.2,15),trace='none')
		dev.off()
		}
}

	
if(args$getSNVs && !args$predictTestSet){	
	allData = allSNPs[which(rownames(allSNPs) %in% typed ),]
	total_classes = classes
	splsda.model <- splsda(t(allData), factor(t(total_classes)), ncomp=1, keepX=optK)
	SNVs = splsda.model$loadings$X
	all_snv_names = rownames(SNVs)
	snvs = matrix(, nrow=optK, ncol=3)
	index = 1
	for(i in 1:length(SNVs)){
		if (abs(SNVs[i]) > 0){
			snvs[index,1:3] = c(all_snv_names[i], SNVs[i], names[all_snv_names[i],2])
			index = index + 1
		}
	}
	snvs[,2] = as.numeric(snvs[,2])
	topSNVs = snvs[order(abs(as.numeric(snvs[,2])), decreasing = T),]
	write(t(topSNVs), ncol=3, file=paste(name, "_top_SNVs.txt", sep=""), sep="\t")
	write(t(topSNVs), ncol=3, file=paste(name, "_top_SNVs.csv", sep=""), sep=",")
	
	if(args$produceFigures){
		snv_names = c()
		afs = c()
		classOutcome = c()
		numTop = min(optK, 15)
		for(e in 1:numTop){
			for(s in 1:ncol(allSNPs)){
				
				snv_name = topSNVs[e,1]
				af = allSNPs[snv_name,s]
				if(classes[s] == 1){
					classOutcome = c(classOutcome, args$nameClass1)	
				}else{
					classOutcome = c(classOutcome, args$nameClass2)
				}
				snv_names = c(snv_names, snv_name) 
				afs = c(afs, af)
				
			}
		}
		
		df = data.frame(snp = factor(snv_names, levels=unique(snv_names)), af = afs, Outcome=classOutcome)
		df = na.omit(df)
		pdf(paste(name,"_af_boxplot.pdf",sep=""), 7,3)
		plot4 = ggplot(df, aes(x=Outcome,y=af)) + geom_boxplot( frame.plot=FALSE,axes=FALSE,outlier.size=0) + facet_grid(. ~ snp) + facet_wrap( ~ snp, ncol=min(ceiling(numTop/2), 8))
		plot4 = plot4 +  geom_point(data=subset(df,af>0),size=1.25, aes(colour=Outcome),position = position_jitter(width = 0.1)) + scale_color_brewer(palette="Set1")
		plot4 = plot4 + theme_bw(base_size=5.5) 
		plot4 = plot4 + scale_y_continuous(breaks=c(0,.5,1))
		plot4 = plot4 + ylab("Allele Fraction")
		plot4 = plot4 + xlab(NULL) + ggtitle(paste("Top ", name, " SNVs",sep=""))
		plot4 = plot4 + theme(axis.title.x=element_blank(),
													axis.text.x=element_blank(),
													axis.ticks.x=element_blank(), legend.position="bottom")
		plot4
		print(plot4)
		dev.off()
		
		pdf(paste(name,"_heatmap.pdf", sep=""), 8,5)
		data = allSNPs[which(rownames(allSNPs) %in% topSNVs[1:numTop,]),]
		SNVs = rownames(data)
		samps = colnames(allSNPs)
		data = matrix(as.numeric(unlist(data)), nrow=nrow(data))
		colnames(data) = samps
		rownames(data) = t(SNVs)
		heatmap.2(data, key.title="NA", key.xlab="Allele Fraction",density.info='none',na.color='black', col = colorpanel(100,"white","blue"),margins=c(4.2,15),trace='none')
		dev.off()
	}	
}



if(args$permutationTests && !args$predictTestSet){

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

	permIter = as.numeric(args$permIter)
	write(paste("Running ",permIter," permuation tests.",sep=""),file=paste(name, ".log", sep=""), append=T)	

	ls <- foreach (iter=1:permIter, .packages=c('mixOmics','pROC')) %dopar% {
	#	for(iter in 1:numIter){
		write(paste("Permutation iteration: ",iter, sep=""), file=paste(name, ".log", sep=""), append=T)
		if(is.null(args$NfoldCV)){
		results = get_predictions(skips, allSNPs[, sample(c(seq(1,ncol(allSNPs),1)))], classes, optK)
		results$AUC[2]
		}else{
			
			allSNPs_permute = allSNPs[, sample(c(seq(1,ncol(allSNPs),1)))]
				#Setting up unstratified NfoldCV for internal CV
			indices = 1:ncol(allSNPs_permute)
			numPer = floor(ncol(allSNPs_permute)/args$NfoldCV)
			
			rem = ncol(allSNPs_permute) %% args$NfoldCV
			internal.skips = matrix(NA, nrow=args$NfoldCV, ncol=numPer+1)
			it = 0
			rando = sample(1:ncol(allSNPs_permute))
			for(nf in 1:args$NfoldCV){
				if(nf <= rem){
					
					internal.skips[((it*args$NfoldCV)+nf),] = rando[((numPer+1)*nf-numPer):(nf*(numPer+1))]
					
				}else{
					internal.skips[((it*args$NfoldCV)+nf),] = c(rando[(numPer*nf-numPer+1+rem):(nf*(numPer)+rem)], NA)
				}
			}
			
			internal.skips = data.frame(internal.skips)
				
			results = get_predictions(internal.skips, allSNPs_permute, classes, optK)
			write(results$AUC[2], ncol=1, file=paste(name, "_permutation_AUCs.txt", sep=""), sep="\t", append=T)
			results$AUC[2]
			
		}
	}

	aucs = unlist(ls)

	p_val = length(which(aucs >= compare_auc))/length(aucs)
	write(paste("Perm p-value: ",p_val, sep=""), file=paste(name, ".log", sep=""), append=T)
	write(p_val, ncol=1, file=paste(name,"_p_val.txt", sep=""),sep="\t")

}
	




