function fitness = custom_fitness(member,ytrain, xtrainSubsetSNP, num_OPLS_fact,desiredNGenes)
xtrain = xtrainSubsetSNP;

keepInxs = find(member(1:end) ==1);
num = length(keepInxs);

m = mean(xtrain(:,keepInxs));  %finds the mean of all "on" genes for DE
mean_Y=mean(ytrain);  %finds mean ytrainSubsetAS which is an array of 1's and 0's

Xres = bsxfun(@minus,xtrain(:,keepInxs), m);%subtract the mean from each "on" gene in xtrainSubsetAS. possibly some sort of residual?
Xres = double(Xres); %4added to avoid error Warning: Matrix is singular to working precision.

Yres = ytrain - mean_Y; %subracts mean of 1's and 0's from all
Yres = double(Yres); %added to avoid error Warning: Matrix is singular to working precision.

predict_opls_custom = @(xtrain,ytrain,xtest) (predict_opls(xtrain,ytrain,xtest,num_OPLS_fact));

cvMse = crossval('mse',Xres,Yres,'predfun',predict_opls_custom);

fitness = []; %creats fitness score type thing
fitness(1) = cvMse; %(1-stats.R2_Y); %some sort of statistic returned
fitness(2) = abs(num-desiredNGenes); % total number of genes that are "on"

