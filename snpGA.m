warning off;
%matlabpool open local % if you want to use parallel computing tools in MATLAB, uncomment this line

% Paths that must be set prior to running this script:
path_to_opls = 'C:/Users/andersonp/Documents/GitHub/opls/matlab'; % Available at https://github.com/Anderson-Lab/OPLS
data_dir = 'C:/Users/andersonp/Desktop/exonic_every_combo_Diffsets/exonic_every_combo_DiffSets/exonic_every_combo_top25_diff'; % This must contain the cross-validation sets
export_dir = '';

addpath(path_to_opls)

% Designates method, data & selection f(n) to use in the GA
selection = @selectiontournament; % Possible choices are: @selectionremainder,@selectionuniform,{@selectionstochunif},@selectionroulette,@selectiontournament
cross = @crossovertwopoint; % Possible choices are: @crossoverheuristic,{@crossoverscattered},{@crossoverintermediate}*
popSize=300;

            % @crossoversinglepoint,@crossovertwopoint,@crossoverarithmetic

% Initializes the "necessary" arrays (note the quotation marks around necessary)
results = {};
numCorrect = {};
numIncorrect = {};
numTimesInTop = {};
avgRank = {};
yPredTrain = {};
yPredTest = {};
yTrain = {};
topN = 100;
trainIndexes = {};
testIndexes = {};
allPops = {};
fvalues = {};
exits = {};
allOuts = {};
allScores = {};

% Need to include frequency of the gene showing up in top 20 (N)
nIterations = 1; % Number of times to run the GA
avgNumVariables = 0;
avgFracASVariables = 0;
desiredNFeatures = 10; % Target number of features
randomness=10; % Randomness specified while creating initial population

gotWrong = [];
allPred = [];
classified = [];
bestFitnessScores = [];
bestFitnessSizes = [];
worstFitnessScores = [];
worstFitnessSizes = [];
Y_actual = [];
Y_real_preds = [];
totalCount = 1;
num_OPLS_fact = 1;  %%% number of orthogonal components in opls model

% These should all come from a metadata file
numClass1 = 11; % This needs to be read in as part of the data file
numClass2 = 9;
ytrainSubsetSNP = [0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1];
ytrainSubsetSNP = ytrainSubsetSNP';
ytestDE = [0,1];
ytestDE = ytestDE';% Iterates through the GA

for z1 = 1:((numClass1))
    for z2 = 1:((numClass2))
        xtrainSubsetSNP = csvread(strcat(data_dir,'/training.data.',num2str(z1-1),'.',num2str(z2-1),'.csv'));
        xtrainSubsetSNP = xtrainSubsetSNP';
        xtestSNP = csvread(strcat(data_dir,'/test.data.',num2str(z1-1),'.',num2str(z2-1),'.csv'));
        xtestSNP = xtestSNP';
        gene_names = dataset('file', strcat(data_dir,'/feature.names.',num2str(z1-1),'.',num2str(z2-1),'.csv'), 'delimiter', ',','ReadVarNames','off');

        
        bestScore = 1000;
        bestSet = [];
        predicted = 0;
        
        leaveOutIndex = [z1,z2+11];
        
        classified = [classified;z1];
        classified = [classified;z2+11];
                                
        % Iterates through the GA
        for n = 1:nIterations            
            % Initializes the fitness funcgtion and population information
            %fitness = @(member) (finalFitness(member,xtrainDE,ytrainDE, xtrainSubsetAS,ytrainSubsetAS, xtrainSNP, xtestSNP, xtestDE, ytestDE, xtestAS, ytestAS ,gene_namesDE,gene_namesAS,a,desiredNFeatures,method, data, xtrainSubsetSD, xtestSD,kernel));
            fitness = @(member) (custom_fitness(member,ytrainSubsetSNP, xtrainSubsetSNP, num_OPLS_fact,desiredNFeatures));
            
            dataSize = size(xtrainSubsetSNP,2);
            initPop = zeros(popSize,dataSize);
            
            % Creates the starting population
            for p = 1:popSize
                ixs = randperm(dataSize);
                initPop(p,ixs(1:round(desiredNFeatures+rand*randomness))) = 1;
            end
            
            % Sets GA options and runs the GA
            options = gaoptimset('PopulationType','bitstring','initialpopulation',initPop, 'SelectionFcn', selection, 'CrossoverFcn', cross, 'populationsize',popSize);
            [x,fval,exitflag,output,population,scores]= gamultiobj(fitness,dataSize,[],[],[],[],[],[],options);
            
            % Stores the results from the GA run            
            geneSubset = [];
            counts = [];            
            for i = 1:popSize                
                set = gene_names(find(population(i,:)),1);
                geneSubset = [geneSubset;set];
                counts(i) = length(find(population(i,:)));
                
            end
            totalGeneSubset = [totalGeneSubset; geneSubset];
            
            [geneFreqNames, geneFreqNums] = count_unique(geneSubset.Var1);
            
            
            bestChrInx = find(scores(:,1) == min(scores(:,1)));
            
            bestChrMin = min(counts(bestChrInx));
            
            best = find(counts(bestChrInx) == bestChrMin);
            best = bestChrInx(best(1));
            bestIndx = find(population(best,:));
            bestGenes = gene_names(find(population(best,:)),1);
            
            bestChrSubset = [bestChrSubset;bestGenes];
            
            
            bestFitnessScores = [bestFitnessScores;scores(best,1)];
            bestFitnessSizes = [bestFitnessSizes;scores(best,2)];
            
            worstChrInx = find(scores(:,1) == max(scores(:,1)));
            
            worstChrMin = max(counts(worstChrInx));
            
            worst = find(counts(worstChrInx) == worstChrMin);
            worst = worstChrInx(worst(1));
            
            worstFitnessScores = [worstFitnessScores;scores(worst,1)];
            worstFitnessSizes = [worstFitnessSizes;scores(worst,2)];
        end
        
        %End GAs for each sample
        
        
        TrainSubset = xtrainSubsetSNP(:,bestIndx);
        TestSubset = xtestSNP(:,bestIndx);
        subset150{totalCount} = bestIndx;
        totalCount = totalCount + 1;
        
        if (strcmp(method,'opls'))
            [model, stats] = opls(TrainSubset, ytrainSubsetSNP, a);
            [t,t_ortho,Y_pred] = apply_opls_model(TrainSubset,ytrainSubsetSNP,model,TestSubset);
            Y_real_preds = [Y_real_preds;Y_pred];
            Y_pred = round(Y_pred);
            for i = 1:length(Y_pred)
                if (Y_pred(i) == ytestDE(i))
                    predicted = predicted + 1;
                else
                    gotWrong = [gotWrong;leaveOutIndex(i)];
                end
            end
        end
        
        
        allPred = [allPred;predicted];
    end
end
save('subset150DE','subset150');
save('ROCSNPWORKSPACE');

allPred = sum(allPred) / length(classified);
[wrongIndx, wrongCounts] = count_unique(gotWrong);
[classifiedNums, classifiedCounts] = count_unique(classified);
[topGeneFreqNames, topGeneFreqNums] = count_unique(bestChrSubset.Var1);

ROC_labels = classified;
for i = 1:length(ROC_labels)
    if (ROC_labels(i) < 12)
        ROC_labels(i) = 0;
    else
        ROC_labels(i) = 1;
    end
end

ROC_scores = Y_real_preds;

df = dataset({},{},'VarNames',{'Sample','Y_Predicted'});
df.Sample = classified;
df.Y_Predicted = Y_real_preds;
export(df,'file','ROC_SNP_Y_Pred_Values.txt');

% TODO: Need to change this and other outputs to use a prefix, so the user
% can easily change where the data should to exported

% EXPORTS ALL GENE FREQUENCIES OUT INTO TXT FILES --> USED PRIMARILY FOR
% CLUSTER RUNS

df = dataset({},{},'VarNames',{'ClassifiedNums','ClassfiedCounts'});
df.ClassifiedNums = classifiedNums;
df.ClassifiedCounts = classifiedCounts;
export(df,'file','ROC_SNP_classified_freq.txt');

df = dataset({},{},'VarNames',{'wrongSamples','wrongCounts'});
df.wrongSamples = wrongIndx;
df.wrongCounts = wrongCounts;
export(df,'file','ROC_SNP_wrong_frequencies.txt');

df = dataset({}, 'VarNames',{'PredictPercent'});
df.PredictPercent = allPred;
export(df,'file','ROC_SNP_predict_percent.txt');

df = dataset({},{},'VarNames',{'TopGeneNames','TopGeneNums'});
df.TopGeneNames = topGeneFreqNames;
df.TopGeneNums = topGeneFreqNums;
export(df,'file','ROC_SNP_top_ALL_frequencies.txt');

bestFitnessAverage = mean(bestFitnessScores);
tdf = dataset({},{},{},'VarNames',{'bestFitnessValues','bestFitnessSizes','bestFitnessAvg'});
tdf.bestFitnessValues = bestFitnessScores;
tdf.bestFitnessSizes = bestFitnessSizes;
tdf.bestFitnessAvg = bestFitnessAverage;
export(tdf,'file','ROC_SNP_best_fitness_values.txt');

worstFitnessAverage = mean(worstFitnessScores);
tdf = dataset({},{},{},'VarNames',{'worstFitnessValues','worstFitnessSizes','worstFitnessAvg'});
tdf.worstFitnessValues = worstFitnessScores;
tdf.worstFitnessSizes = worstFitnessSizes;
tdf.worstFitnessAvg = worstFitnessAverage;
export(tdf,'file','ROC_SNP_worst_fitness_values.txt');

[ROC_X,ROC_Y,ROC_T,ROC_AUC] = perfcurve(ROC_labels,ROC_scores,1,'nboot',1000);

df = dataset({},{},{},'VarNames',{'ROC_X','ROC_Y','ROC_AUC'});
df.ROC_X = ROC_X;
df.ROC_Y = ROC_Y;
df.ROC_AUC = ROC_AUC;
export(df,'file','ROC_SNP_FigStats.txt');

f = figure;

plot(ROC_X,ROC_Y);
saveas(f,'ROC_SNP_Plot.fig');

matlabpool close



