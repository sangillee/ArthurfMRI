% fmriprepBetaSeries
% function to extract single trial betas from fmriprep data
% 2021-02-06: created
% 2021-02-16: introduced option struct to specify regressions

function [coef,mask,outcome] = fmriprepBetaSeries(fmriprepdir,subjname,runnum,task,session,trialREG,TR,method,nuisreg,opts)
% if subject name already contains sub- prefix, remove it
if strcmp(subjname(1:4),'sub-')
    subjname = subjname(5:end);
end

% load stuff
[img,mask,var] = helper_fmriprep_loadVars(fmriprepdir,subjname,runnum,task,session);
covariates = helper_fmriprep_defaultCovariates(var); % load covariates

% create regressors
nvol = size(img,1); X = nan(nvol,length(trialREG)); nT = size(X,2);
for i = 1:nT
    X(:,i) = simBOLD(TR,nvol,trialREG{i});
end

% create nuisance regressors
if nargin > 8
    nuisX = nan(nvol,length(nuisreg));
    for i = 1:size(nuisX,2)
        nuisX(:,i) = simBOLD(TR,nvol,nuisreg{i});
    end
else
    nuisX = nan(nvol,0);
end

% remove first volume since first volume does not have motion regressor from previous volume
X(1,:) = []; img(1,:) = []; covariates = helper_fmriprep_makeFullRank(covariates(2:end,:)); nuisX(1,:) = [];
X = X - mean(X); covariates = covariates - mean(covariates); img = img - mean(img); nuisX = nuisX - mean(nuisX); % mean centering to remove intercept

if strcmp(method,'LSA')
    coef = helper_fmriprep_massRegression(double([X,covariates,nuisX]),double(img));
    coef = coef(1:nT,:);
elseif strcmp(method,'LSS')
    coef = nan(nT,size(img,2));
    Xsum = sum(X,2);
    for i = 1:nT
        Xint = X(:,i); Xrest = Xsum-Xint;
        tempcoef = helper_fmriprep_massRegression(double([Xint,Xrest,covariates,nuisX]),double(img));
        coef(i,:) = tempcoef(1,:);
    end
else
    error('unknown beta regression type')
end

end