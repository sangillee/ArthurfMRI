% searchlight
% Arthur Lee
% function that performs searchlight analysis
% Usage: measure = searchlight(betamaps,CVfold,searchlightInd)
% 2021-01-31 Written
% 2021-02-06 Included AUC in metric

function measure = searchlight(betamaps,Y,CVfold,searchlightInd,metric)
% betamaps = n by v matrix of beta coefficients
% CVfold: a binary matrix where each column marks the testing fold as 1 and training fold as 0
if size(CVfold,2) == 1
    temp = CVfold; uniqfold = unique(temp);
    CVfold = zeros(length(Y),length(uniqfold));
    for i = 1:length(uniqfold)
        CVfold(temp==uniqfold(i),i) = 1;
    end
end

[n,v] = size(betamaps);
measure = nan(size(CVfold,2),v);
for i = 1:v % iterating through each voxel
    i/v
    searchind = searchlightInd(:,i); searchind = searchind(~isnan(searchind));
    X = betamaps(:,searchind); remind = (std(X)==0) | isnan(std(X)) ; X(:,remind) = [];
    if ~isempty(X)
        for j = 1:size(CVfold,2) % iterating through each fold
            trainind = CVfold(:,j)==0;
            trainX = X(trainind,:); trainY = Y(trainind);
            testX = X(~trainind,:); testY = Y(~trainind);
            measure(j,i) = dopred(trainX,trainY,testX,testY,metric);
        end
    end
end
end

function perf = dopred(trainX,trainY,testX,testY,metric)
B = regress(trainY,[ones(size(trainY)),trainX]);
pred = [ones(size(testY)),testX]*B;
if strcmp(metric,'Pearson') || strcmp(metric,'Spearman')
    if std(pred) == 0
        perf = 0;
    else
        perf = corr(pred,testY,'type',metric);
    end
elseif strcmp(metric,'AUC')
    perf = localAUC(binarize(testY)==1,pred);
else
    error('unknown metric')
end
end

function auc = localAUC(labels,pred)
assert(size(labels,1)==size(pred,1),'number of rows do not match for inputs')
assert( islogical(labels) ,'labels input should be logical');
assert( size(labels,2)==1, 'labels should not be a matrix');
n = size(labels,1); num_pos = sum(labels); num_neg = n - num_pos;
auc = nan(1,size(pred,2));
if (num_pos>0 && num_pos < n)
    ranks = tiedrank(pred);
    auc = ( sum( ranks(labels,:) ) - num_pos * (num_pos+1)/2) / ( num_pos * num_neg);
end
end

function vect = binarize(vect)
% convert a vector that has two values into fals and true
vals = unique(vect);
if length(vals) ~= 2
    error('more than two categories detected')
end
cat1 = vect==vals(1); cat2 = vect==vals(2);
vect(cat1) = 0; vect(cat2) = 1;
end