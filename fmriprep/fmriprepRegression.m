% fmriprepRegression
% 2021-02-06: added a_comp_cor regressors into pipeline. separated covariate loading into a function
% 2021-02-16: introduced option struct to specify regressions

function [coef,mask,outcome] = fmriprepRegression(fmriprepdir,subjname,runnum,task,session,REG,TR,opts)
if nargin< 8
    opts = helper_fmriprep_regoptions; % loading default options
end

% load stuff
[img,mask,var] = helper_fmriprep_loadVars(fmriprepdir,subjname,runnum,task,session);
covariates = helper_fmriprep_covariates(var,opts.covar); % load covariates

% smoothing, if requested
if opts.FWHM > 0
    img = smoothfMRI(img,opts.FWHM,2,mask);
end

[coef,outcome] = fitfMRI(REG,double(img),covariates,TR,opts.HP,opts.measure,opts.HRF);
end