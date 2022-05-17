% fmriprepRegression
% 2021-02-06: added a_comp_cor regressors into pipeline. separated covariate loading into a function
% 2021-02-16: introduced option struct to specify regressions
% 2022-03-27: removed default options. The user must supply a mask to work with
% 2022-04-22: included ability to run many regressions with one image depending on the shape of the REG structure

function [coef,stats] = fmriprepRegression(fmriprepdir,subjname,runnum,task,session,REG,TR,opts,mask)
% check option input for errors
helper_fmriprep_regoptions(opts);

% load stuff
[img,var] = helper_fmriprep_loadVars(fmriprepdir,subjname,runnum,task,session);
img = double(img);
covariates = helper_fmriprep_covariates(var,opts.covar); % check option structure and load covariates

% smoothing, if requested
if opts.FWHM > 0
    img = smoothfMRI(img,opts.FWHM,2);
end

img = reshape(img,size(img,1)*size(img,2)*size(img,3),size(img,4))'; % re-shaping the dimensions
img = img(:,mask(:)==1);

nrep = size(REG,3);
if nrep == 1
    [coef,stats] = fitfMRI(REG,img,covariates,TR,opts.HP,opts.measure,opts.HRF,opts.drop);
else
    coef = cell(1,nrep); stats = cell(1,nrep);
    for i = 1:nrep
        [coef{i},stats{i}] = fitfMRI(REG(:,:,i),img,covariates,TR,opts.HP,opts.measure,opts.HRF,opts.drop);
    end
end
end