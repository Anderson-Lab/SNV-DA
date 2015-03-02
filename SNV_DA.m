%warning off;
%matlabpool open local % if you want to use parallel computing tools in MATLAB, uncomment this line
function snv_run = SNV_DA(data_dir, export_dir, popSize,numClass1, numSamps, desiredN, randAmount)
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
desiredNFeatures = desiredN; % Target number of features
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


numClass2 = numSamps - numClass1;

ytrainSubset = [zeros(numClass1-1,1);ones(numClass2-1,1)];
ytest = [0,1];
ytest = ytest';% Iterates through the GA
GAcount = 1;
first = 1;
totalGeneSubset = [];
bestChrSubset = [];
bestFitnessScores = [];
numModels = 0;
for z1 = 1:((numClass1))
    for z2 = 1:((numClass2))
        
        present = exist(strcat(data_dir,'/training.data.',num2str(z1-1),'.',num2str(z2-1),'.csv'), 'file');
        if present ~= 2
            continue
        end
        numModels = numModels + 1;
        disp(strcat('GA: ',num2str(GAcount))); 
        GAcount = GAcount + 1;
        xtrainSubset = csvread(strcat(data_dir,'/training.data.',num2str(z1-1),'.',num2str(z2-1),'.csv'));
        xtrainSubset = xtrainSubset';
        xtest = csvread(strcat(data_dir,'/test.data.',num2str(z1-1),'.',num2str(z2-1),'.csv'));
        xtest = xtest';
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
           
            fitness = @(member) (custom_fitness(member,ytrainSubset, xtrainSubset, num_OPLS_fact,desiredNFeatures));
            
            dataSize = size(xtrainSubset,2);
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
        
        
        TrainSubset = xtrainSubset(:,bestIndx);
        TestSubset = xtest(:,bestIndx);
        subset150{totalCount} = bestIndx;
        totalCount = totalCount + 1;
        
	
        [model, stats] = opls(TrainSubset, ytrainSubset, num_OPLS_fact);
        [t,t_ortho,Y_pred] = apply_opls_model(TrainSubset,ytrainSubset,model,TestSubset);
         Y_real_preds = [Y_real_preds;Y_pred];
        Y_pred = round(Y_pred);
        for i = 1:length(Y_pred)
        if (Y_pred(i) == ytest(i))
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

group1_wrong = 0;
group2_wrong = 0;
for i = 1:length(wrongIndx)
    if (wrongIndx(i) <= numClass1)
        group1_wrong = group1_wrong + wrongCounts(i);
    else
        group2_wrong = group2_wrong + wrongCounts(i);
    end
end

group1_sens = (numModels - group1_wrong) / numModels;
group2_sens = (numModels - group2_wrong) / numModels;

df = dataset({},{},'VarNames',{'Top_SNVs','Model_Frequenxt'});
df.TopGeneNames = topGeneFreqNames;
df.TopGeneNums = topGeneFreqNums/numModels;
export(sortrows(df,'TopGeneNums','descend'),'file',strcat(export_dir,'/top_SNVs.txt'));

ROC_scores = Y_real_preds;
[ROC_X,ROC_Y,~,ROC_AUC] = perfcurve(ROC_labels,ROC_scores,1,'nboot',1000);

f = figure;

[ROC_X,ROC_Y] = perfcurve(ROC_labels,ROC_scores,1);
plot(ROC_X,ROC_Y);
saveas(f,strcat(export_dir,'/ROC_plot.fig'));

report=fopen(strcat(export_dir,'/model_stats_report.txt'),'w');
fprintf(report, '%6s\t%12s\t%6s\t%12s\n','Prediction_Accuracy', 'AUC[95% CI]', 'Group1_Sensitivity', 'Group2_Sensitivity');
fprintf(report, '%.4f\t%.4f[%.4f-%.4f]\t%.4f\t%.4f\n', allPred, ROC_AUC, group1_sens, group2_sens);
fclose(report);






