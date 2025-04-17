% helper_fmriprep_loadVars
% function to take care of loading images, masks, and covariates
% 2021-02-06: created
% 2021-09-02: incorporated missing runnum
% 2022-03-22: moved smoothing in here instead of fmriprepRegression to smooth before applying the mask instead of after
%             added dilation to run-level mask to avoid missing data at the group-level
% 2022-03-27: removed using run-specific brain mask. moved smoothing back outside the shell
% 2022-05-04: included grand mean scaling
% 2023-04-02: more robust file loading through getonlyfile function
% 2025-04-15: Cleaning up for packaging

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
img = Niftiopen(getonlyfile([address,'_space*desc-preproc_bold.nii*'])); % img = load_nii(getonlyfile([address,'_space*desc-preproc_bold.nii*'])); img = img.img;

% grand mean scaling
mask = Niftiopen(getonlyfile([address,'_space*desc-brain_mask.nii*'])); % mask = load_nii(getonlyfile([address,'_space*desc-brain_mask.nii*'])); mask = mask.img;
tempimg = reshape(img,[],size(img,4)); tempimg = tempimg(mask(:)==1,:);
grandmean = abs(median(tempimg(:)));
img = img .* (10000/grandmean);

var = readtable(getonlyfile([address,'_desc-confounds*.tsv']),'FileType','text','TreatAsEmpty','n/a');
end

function filename = getonlyfile(filename)
    % when you know the approximate file name but only want to load the file if there's only one instance of it
    out = dir(filename);
    if length(out) ~= 1
        error(['Error in loading file ',filename,'. Multiple files detected.'])
    end
    filename = [out.folder,'\',out.name];
end