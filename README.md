# SNV-DA
 
Rscript snvDA.R -h
Usage: snvDA.R [options]

SNV-DA is used to create and evaluate sPLS-DA models to identify single nucleotide variations (SNVs) that accurately classify phenotype. The steps of the pipeline include:
1) Finding optimal value of number of features to be selected in the model (K)
2) Evaluating the model using cross-validations and optimal value of K
3) Find and rank the K features that are correlated with predictive accuracy
4) Permutation tests to determine if model is discriminate towards the true grouping of labels.


The script has several parameters that allows the user to design their own analysis:
1) The number of non-zero values by which to filter the matrix
2) The number of NA values allowed for each SNV
3) The type of SNV to include in the model
4) The range and/or values of K tested.
4) Cross-validation design (number of samples to take out from each group, the stratification of test samples [force test samples to be pulled equally from each group], and the number of cross-validations).

Happy hunting!

Options:
	-O, --findOptimalK
		If flagged, SNV-DA will run nest cross-validations to determine optimal K.

	-J, --evaluatePerformance
		If flagged, SNV-DA will evaluate performance of model the value of K found from running --findOptimalK or user supplied --optimalK.

	-P, --permutationTests
		If flagged, SNV-DA will run permutation tests by randomly permuting sample classes then running cross-validations. Tests are run --permIter number of times. Performances are then compared to the value of AUC found by --evaluatePerformance of user supplied --AreaUnderCurve.

	-x, --predictTestSet
		If flagged, SNV-DA will train a model with optimal K on all samples in the main matrix and then predict the labels of the second matrix.

	-C, --getSNVs
		If flagged, SNV-DA will ONLY run sPLS-DA on all samples to identify K top predictive features

	-W NFOLDCV, --NfoldCV=NFOLDCV
		If specified, N-fold cross-validations will be performed by partitioning the input matrix into NfoldCV setsi during the optimization of K and performance evaluation.

	-I ITERNFOLDCV, --iterNfoldCV=ITERNFOLDCV
		If NfoldCV is specified, NfoldCV will be run iterNFoldCV times during performance evaluation  and finding of optimal K (default=1). 

	-V INTERNALNFOLDCV, --internalNfoldCV=INTERNALNFOLDCV
		If specified, N-fold cross-validations will be performed during the interval CV while finding optimal K by partitioning the input matrix into iternalNfoldCV setsi during the optimization of K and performance evaluation. If left unspecified, leave-one-out crossvalidations will be used for the interval cross-validation.

	-j INTERNALITERNFOLDCV, --internalIterNfoldCV=INTERNALITERNFOLDCV
		If internalNfoldCV is specified, internalNfoldCV will be run internalIterNFoldCV times during performance evaluation and finding of optimal K (default=1).

	-M MATRIX, --matrix=MATRIX
		SNVM input (csv)

	-y TESTMATRIX, --testMatrix=TESTMATRIX
		matrix to predict labels (csv)

	-D SIZECLASS1, --sizeClass1=SIZECLASS1
		The number of samples of class 1. (Used to find the number of samples of both classes).

	-L MINK, --minK=MINK
		The minimum value of K tested.

	-B MAXK, --maxK=MAXK
		The maximum value of K tested.

	-E EVERYKTH, --everyKth=EVERYKTH
		Specifies every Nth K tested.
Or use:

	-R RANGE, --range=RANGE
		A comma separated list of K's to be added to the Ks tested, e.g. 5,10,15,20,25,50,75,100,250,300.

	-Z NUMNONZEROS, --numNonZeros=NUMNONZEROS
		Limit SNVs to those that have at least X number of samples with allele fractions greater than zero (default = 3). Decreases runtime. Suggested value is a 1/3 the number of samples in the smaller group.

	-N NUMNAS, --numNAs=NUMNAS
		Limit SNVs to those that have at most X amount of samples with 'NA' values (default=0).

	-G INCLUDE, --include=INCLUDE
		SNVs with this substring are included, e.g. "exonic" for nonsynonymous and synonymous exonic SNVs

	-X EXCLUDE, --exclude=EXCLUDE
		SNVs with this substring are excluded, e.g. ncRNA to remove intronic_ncRNA

	-A STUDYNAME, --studyName=STUDYNAME
		One word name of the analysis used in output file names and figures.

	-U NAMECLASS1, --nameClass1=NAMECLASS1
		Name of class 1 (default = Class_1).

	-H NAMECLASS2, --nameClass2=NAMECLASS2
		Name of class 2 (default = Class_2).

	-T THREADS, --threads=THREADS
		Number of threads for cross-validations (default = 1).

	-K OPTIMALK, --optimalK=OPTIMALK
		Used when --evaluatePerformance or --permuationTests is flagged. Specifies the value of K used for cross-validation (default = 500).

	-q PERMITER, --permIter=PERMITER
		Used when --permuationTests is flagged. Specifies the number of permutation tests to be run

	-a AREAUNDERCURVE, --AreaUnderCurve=AREAUNDERCURVE
		Performance of true model. Compared to permuted models when --permutationTests are run. If unflagged, value is determined by --evaluatePerformance

	-F, --produceFigures
		If flagged, allele fraction box plots, heatmaps, and kernal density figures are produced.

	-h, --help
		Show this help message and exit
