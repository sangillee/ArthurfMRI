% helper_fmriprep_loadVars
% function to take care of loading images, masks, and covariates
% 2021-02-06: created
% 2021-09-02: incorporated missing runnum
% 2022-03-22: moved smoothing in here instead of fmriprepRegression to smooth before applying the mask instead of after
%             added dilation to run-level mask to avoid missing data at the group-level
% 2022-03-27: removed using run-specific brain mask. moved smoothing back outside the shell

function [img,var] = helper_fmriprep_loadVars(fmriprepdir,subjname,runnum,task,session)
% if subject name already contains sub- prefix, remove it
if length(subjname)>3 && strcmp(subjname(1:4),'sub-')
    subjname = subjname(5:end);
end

% set up fmriprep file directory
if isempty(session)
    address = [fmriprepdir,'sub-',subjname,'/func/sub-',subjname,'_task-',task];
else
    address = [fmriprepdir,'sub-',subjname,'/ses-scan',num2str(session),'/func/sub-',subjname,'_ses-scan',num2str(session),'_task-',task];
end
if ~isempty(runnum)
    address = [address,'_run-',num2str(runnum)];
end

% load images
img = load_nii([address,'_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz']);
img = img.img;

if isfile([address,'_desc-confounds_regressors.tsv'])
    var = readtable([address,'_desc-confounds_regressors.tsv'],'FileType','text','TreatAsEmpty','n/a'); % older version of fmriprep
else
    var = readtable([address,'_desc-confounds_timeseries.tsv'],'FileType','text','TreatAsEmpty','n/a');
end
end