% fmriprepBetaSeries
% function to extract single trial betas from fmriprep data
% trialREG is a 3 column matrix where each row corresponds to each trial
% 2021-02-06: created
% 2021-02-16: introduced option struct to specify regressions

function [coef,mask,outcome] = fmriprepBetaSeries(fmriprepdir,subjname,runnum,task,session,trialREG,TR,method,nuisreg,opts)
if nargin < 10
    opts = helper_fmriprep_regoptions; % loading default options
end
if nargout == 3 && strcmp(opts.measure,'nothing')
    error('must specify measure in option struct')
end

% load stuff
[img,mask,var] = helper_fmriprep_loadVars(fmriprepdir,subjname,runnum,task,session);
covariates = helper_fmriprep_covariates(var,opts.covar); % load covariates

% smoothing, if requested
if opts.FWHM > 0
    img = smoothfMRI(img,opts.FWHM,2,mask);
end

nT = size(trialREG,1); v = size(img,2);
if strcmp(method,'LSA')
    REG = cell(1,nT);
    for i = 1:nT
        REG{i} = trialREG(i,:);
    end
    [coef,outcome] = fitfMRI([REG,nuisreg],double(img),covariates,TR,opts.HP,opts.measure,opts.HRF);
    coef = coef(1:nT,:);
    if ~strcmp(opts.measure,'nothing')
        outcome = outcome(1:nT,:);
    end
elseif strcmp(method,'LSS')
    coef = nan(nT,v); outcome = nan(nT,v);
    for i = 1:nT
        REGint = {trialREG(i,:)}; REGnuis = {trialREG([1:(i-1),(i+1):nT],:)};
        [tempcoef,tempoutcome] = fitfMRI([REGint,REGnuis,nuisreg],double(img),covariates,TR,opts.HP,opts.measure,opts.HRF);
        coef(i,:) = tempcoef(1,:);
        if ~strcmp(opts.measure,'nothing')
            outcome(i,:) = tempoutcome(1,:);
        end
    end
else
    error('unknown beta regression type')
end
end