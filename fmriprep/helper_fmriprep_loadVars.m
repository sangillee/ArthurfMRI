% helper_fmriprep_loadVars
% function to take care of loading images, masks, and covariates
% 2021-02-06: created

function [img,mask,var] = helper_fmriprep_loadVars(fmriprepdir,subjname,runnum,task,session)
% if subject name already contains sub- prefix, remove it
if length(subjname)>3 && strcmp(subjname(1:4),'sub-')
    subjname = subjname(5:end);
end

% set up fmriprep file directory
if isempty(session)
    address = [fmriprepdir,'sub-',subjname,'/func/sub-',subjname,'_task-',task,'_run-',num2str(runnum)];
else
    address = [fmriprepdir,'sub-',subjname,'/ses-scan',num2str(session),'/func/sub-',subjname,'_ses-scan',num2str(session),'_task-',task,'_run-',num2str(runnum)];
end

% load images
img = load_nii([address,'_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz']);
img = reshape(img.img,size(img.img,1)*size(img.img,2)*size(img.img,3),size(img.img,4))'; % re-shaping the dimensions
mask = load_nii([address,'_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz']); mask = mask.img; img = img(:,mask(:)==1);

if isfile([address,'_desc-confounds_regressors.tsv'])
    var = readtable([address,'_desc-confounds_regressors.tsv'],'FileType','text','TreatAsEmpty','n/a'); % older version of fmriprep
else
    var = readtable([address,'_desc-confounds_timeseries.tsv'],'FileType','text','TreatAsEmpty','n/a');
end
end