% searchlight function by Arthur Lee
% function that performs searchlight analysis
%
% 2021-01-31 Written
% 2021-02-06 Included AUC in metric
% 2022-07-28 Includes searchlight neighbor function. Defaults to PLS1
% 2023-04-14 Added OLS
% 2024-08-16 added ridge regression

function measure = searchlight(X,Y,CVfold,mask,dist,metric,method)
    if nargin < 7; method = 'PLS1'; end
    if nargin < 6; metric = 'Pearson'; end
    if nargin < 5; dist = 2; end
    assert(ismember(method,{'PLS1','OLS','ridge'}),'unknown method')
    assert(ismember(metric,{'Pearson','AUC'}),'unknown metric')
    assert(size(CVfold,2) == 1,'CV fold should be a vector input')
    assert(all(~isnan(X(:))) && all(~isnan(Y(:))),'there should be no NaN')
    [CVfold,nfold] = generateCVfold(CVfold); % vet CVfold
    neighborList = searchlight_getNeighbors(mask,dist); % get searchlight neighbor list

    measure = nan(nfold,size(X,2));
    for i = 1:1 % iterate through folds
        i
        trainind = CVfold~=i;
        score = doPred(X(trainind,:), Y(trainind), X(~trainind,:),neighborList,method,CVfold(trainind),metric);
        measure(i,:) = getPerf(Y(~trainind),score,metric);
    end
end

function score = doPred(trainX,trainY,testX,neighborList,method,trainCVfold,metric)
    switch method
        case 'PLS1' % using only the first component of Partial Least Squares. This greatly simplifies the calculations and can avoid loops
            CovMap = (trainX-mean(trainX))'*(trainY-mean(trainY)); % this is the model essentially
            predictor = testX .* repmat(CovMap',size(testX,1),1);
            score = cell2mat(cellfun(@(x) sum(predictor(:,x),2), neighborList, 'UniformOutput', false)');
        case 'OLS' % no tuning required
            score = nan(size(testX)); n = size(testX,1);
            for i = 1:size(trainX,2) % iterating through each voxel
                score(:,i) = [ones(n,1),testX(:,neighborList{i})]* regress(trainY,[ones(size(trainY)),trainX(:,neighborList{i})]);
            end
        case 'ridge'
            score = nan(size(testX)); n = size(testX,1);
            for i = 1:size(trainX,2) % iterating through each voxel
                score(:,i) = [ones(n,1),testX(:,neighborList{i})]* localCVRidge(trainX(:,neighborList{i}),trainY,trainCVfold,metric); 
            end
    end
end

function b = localCVRidge(X,Y,CVfold,metric)
    k = [logspace(-5,5,21)];
    [CVfold,nfold] = generateCVfold(CVfold);
    measure = nan(nfold,length(k));
    for i = 1:nfold
        b = localRidge(Y(CVfold~=i),X(CVfold~=i,:),k);
        score = b(1,:) + X(CVfold==i,:)*b(2:end,:);
        measure(i,:) = getPerf(Y(CVfold==i),score,metric);
    end
    [~,ind] = max(nanmean(measure));
    b = localRidge(Y,X,k(ind));
end

function b = localRidge(y,X,k)
    [n,p] = size(X); nk = numel(k);

    % Normalize the columns of X to mean zero, and standard deviation one.
    mx = mean(X); stdx = std(X,0,1);
    idx = find(abs(stdx) < sqrt(eps(class(stdx))));
    if any(idx); stdx(idx) = 1; end
    MX = mx(ones(n,1),:); STDX = stdx(ones(n,1),:);
    Z = (X - MX) ./ STDX;
    if any(idx); Z(:,idx) = 1; end

    % adding pseudo observations having y=0 and X'X = k*I.
    Zplus  = [Z;sqrt(k(1)) * eye(p)]; 
    yplus  = [y;zeros(p,1)];

    % Compute the coefficient estimates
    b = Zplus\yplus;
    if nk>1
        b(end,nk) = 0;
        for j=2:nk
            Zplus(end-p+1:end,:) = sqrt(k(j)) * eye(p);
            b(:,j) = Zplus\yplus;
        end
    end
    b = b ./ repmat(stdx',1,nk); % Put on original scale
    b = [mean(y)-mx*b; b];
end

function perf = getPerf(testY,pred,metric)
    if strcmp(metric,'Pearson')
        perf = corr(pred,testY);
    else
        perf = AUC(testY>0,pred);
    end
end

function [CVfold,nfold] = generateCVfold(inCVfold) % convert CVfold into unique elements from 1~nfold
    uniqfold = unique(inCVfold); nfold = length(uniqfold);
    CVfold = zeros(size(inCVfold));
    for i = 1:nfold
        CVfold(inCVfold==uniqfold(i)) = i;
    end
end