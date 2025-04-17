% fmriprepRegression
% function for running glm on fmriprep-based directory structure
% process:
%   1: input options are checked
%   2: helper_fmriprep_loadVars 1) loads the task brain epi time series
%                               2) loads corresponding brain mask and performs grand mean scaling, and
%                               3) loads all covariates calculated by fmriprep
%   3: image is cast into double precision, and helper_fmriprep_covariates loads covariates selected by the user via opts
%   4: image is smoothed, if requested via opts
%   5: regression is performed via fitfMRI
%
% inputs
%   fmriprepdir: the root directory of fmriprep preprocessed files
%   subjname: subject name, as indicated in the fmriprep folders
%   runnum: if available, the run number of the task. if unavailable, use empty array []
%   task: task name
%   session: if available, the session number. if unavailable, use empty array []
%   REG: 1 x p cell array where each cell contains a three-column FSL-style matrix (first column: onset, second column: duration, third column: modulation)
%        Each cell will create a expected HRF signal and be entered as regressors.
%        See fitfMRI as to how this is used
%   TR: repetition time, in seconds
%   opts: regression options. Check helper_fmriprep_covariates for full array of options.
%   mask: mask to use before regression.
% 
% for outputs, see fitfMRI

% 2021-02-06: added a_comp_cor regressors into pipeline. separated covariate loading into a function
% 2021-02-16: introduced option struct to specify regressions
% 2022-03-27: removed default options. The user must supply a mask to work with
% 2022-04-22: included ability to run many regressions with one image depending on the shape of the REG structure
% 2025-04-15: cleaning up for packaging

function [coef,stats] = fmriprepRegression(fmriprepdir,subjname,runnum,task,session,REG,TR,opts,mask)
% check option input for errors
helper_fmriprep_regoptions(opts);

% load stuff
[img,var] = helper_fmriprep_loadVars(fmriprepdir,subjname,runnum,task,session);
img = double(img);
covariates = helper_fmriprep_covariates(var,opts.covar); % check option structure and load covariates

% smoothing, if requested
if opts.FWHM > 0
    for i = 1:size(img,4)
        img(:,:,:,i) = imgaussfilt3(img(:,:,:,i),opts.FWHM/2.355/2); % convert FWHM to std in voxel units
    end
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