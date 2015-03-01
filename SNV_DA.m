%warning off;
%matlabpool open local % if you want to use parallel computing tools in MATLAB, uncomment this line
function snv_run = SNV_DA(data_dir, export_dir, popSize,numClass1, numSamps, randAmount)
disp('Running SNV-DA')
% Paths that must be set prior to running this script:
path_to_opls = '/Users/DRMRPAUL/repos/OPLS/matlab/'; % Available at https://github.com/Anderson-Lab/OPLS
mkdir(export_dir)
addpath(path_to_opls)
snv_run = 1;
% Designates method, data & selection f(n) to use in the GA
selection = @selectiontournament; % Possible choices are: @selectionremainder,@selectionuniform,{@selectionstochunif},@selectionroulette,@selectiontournament
cross = @crossovertwopoint; % Possible choices are: @crossoverheuristic,{@crossoverscattered},{@crossoverintermediate}*
%popSize=4;

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
randomness=randAmount; % Randomness specified while creating initial population

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
% This needs to be read in as part of the data file
%if first == 1
%    [m,n] = size(xtrainSubsetSNP);
%    numClass2 = numClass1 - m;
%    disp(numClass2);
%end

numClass2 = numSamps - numClass1;

ytrainSubsetSNP = [zeros(numClass1-1,1);ones(numClass2-1,1)];
ytestDE = [0,1];
ytestDE = ytestDE';% Iterates through the GA
GAcount = 1;
first = 1;
totalGeneSubset = [];
bestChrSubset = [];
bestFitnessScores = [];
for z1 = 1:((numClass1))
    for z2 = 1:((numClass2))
        
        present = exist(strcat(data_dir,'/training.data.',num2str(z1-1),'.',num2str(z2-1),'.csv'), 'file');
        if present ~= 2
            continue
        end
            
        disp(strcat('GA: ',num2str(GAcount))); 
        GAcount = GAcount + 1;
        xtrainSubsetSNP = csvread(strcat(data_dir,'/training.data.',num2str(z1-1),'.',num2str(z2-1),'.csv'));
        xtrainSubsetSNP = xtrainSubsetSNP';
        xtestSNP = csvread(strcat(data_dir,'/test.data.',num2str(z1-1),'.',num2str(z2-1),'.csv'));
        xtestSNP = xtestSNP';
        gene_names = dataset('file', strcat(data_dir,'/feature.names.',num2str(z1-1),'.',num2str(z2-1),'.csv'), 'delimiter', ',','ReadVarNames','off');

        
        bestScore = 1000;
        bestSet = [];
        predicted = 0;
        
        leaveOutIndex = [z1,z2+numClass1];
        
        classified = [classified;z1];
        classified = [classified;z2+numClass1];
                                
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
        
	
        [model, stats] = opls(TrainSubset, ytrainSubsetSNP, num_OPLS_fact);
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
        
        
        allPred = [allPred;predicted];
    end
end
save(strcat(export_dir,'/top150'),'subset150');
save(strcat(export_dir,'/SNV_DA_workspace'));

allPred = sum(allPred) / length(classified);
[wrongIndx, wrongCounts] = count_unique(gotWrong);
[classifiedNums, classifiedCounts] = count_unique(classified);
[topGeneFreqNames, topGeneFreqNums] = count_unique(bestChrSubset.Var1);

ROC_labels = classified;
for i = 1:length(ROC_labels)
    if (ROC_labels(i) <= numClass1)
        ROC_labels(i) = 0;
    else
        ROC_labels(i) = 1;
    end
end

ROC_scores = Y_real_preds;

df = dataset({},{},'VarNames',{'Sample','wrongCounts'});
df.Sample = wrongIndx;
df.wrongCounts = wrongCounts;
export(df,'file',strcat(export_dir,'/ROC_SNP_wrong_frequencies.txt'));

df = dataset({}, 'VarNames',{'PredictPercent'});
df.PredictPercent = allPred;
export(df,'file',strcat(export_dir,'/ROC_SNP_predict_percent.txt'));

df = dataset({},{},'VarNames',{'TopGeneNames','TopGeneNums'});
df.TopGeneNames = topGeneFreqNames;
df.TopGeneNums = topGeneFreqNums;
df
export(sortrows(df,'TopGeneNums','descend'),'file',strcat(export_dir,'/ROC_SNP_top_ALL_frequencies.txt'));

[ROC_X,ROC_Y,~,ROC_AUC] = perfcurve(ROC_labels,ROC_scores,1,'nboot',1000);

df = dataset({},{},{},'VarNames',{'ROC_X','ROC_Y','ROC_AUC'});
df.ROC_X = ROC_X;
df.ROC_Y = ROC_Y;
df.ROC_AUC = ROC_AUC;
export(df,'file',strcat(export_dir,'/ROC_SNP_FigStats.txt'));

f = figure;

[ROC_X,ROC_Y] = perfcurve(ROC_labels,ROC_scores,1);
plot(ROC_X,ROC_Y);
saveas(f,strcat(export_dir,'/ROC_SNP_Plot.fig'));

matlabpool close



