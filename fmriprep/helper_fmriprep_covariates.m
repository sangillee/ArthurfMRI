% helper_fmriprep_covariates
% function to take care of loading images, masks, and covariates
% 2021-02-06: created
% 2021-02-16: renamed to helper_fmriprep_Covariates from helper_fmriprep_defaultCovariates
%             now accepts option structure to identify which covariates to include
% 2022-03-27: removed default optioning
%             removed first volume index regressor for t-1 motion regressors

function covariates = helper_fmriprep_covariates(var,opts)
covariates = []; % matrix for covariates

% cosine components
for i = 1:opts.cosines
    covariates = [covariates,loadcovar({['cosine',sprintf('%02d',i-1)]},var)];
end

% a_comp_cor components
for i = 1:opts.a_comp_cor
    covariates = [covariates,loadcovar({['a_comp_cor_',sprintf('%02d',i-1)]},var)]; % adding top CSF+WM components.
end

% csf average
if opts.csf==1
    covariates = [covariates,loadcovar({'csf'},var)];
end

% white matter average
if opts.white_matter ==1
    covariates = [covariates,loadcovar({'white_matter'},var)];
end

% motion params
switch opts.motionparams
    case 0
        motion = [];
    case 6
        motion = loadcovar({'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'},var);
    case 12
        motion = loadcovar({'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z','trans_x_power2','trans_y_power2','trans_z_power2','rot_x_power2','rot_y_power2','rot_z_power2'},var);
    case 24
        motion = loadcovar({'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z','trans_x_power2','trans_y_power2','trans_z_power2','rot_x_power2','rot_y_power2','rot_z_power2'},var);
        motion = [motion,[mean(motion);motion(1:end-1,:)]]; % t-1 motion with mean-imputed first volume
end
covariates = [covariates,motion];

% make full rank
covariates = makeFullRank(covariates);
end

function X = makeFullRank(X)
% remove collinear columns
[n,ncolX] = size(X); [~,R,perm] = qr(X,0);
if isempty(R)
    p = 0;
elseif isvector(R)
    p = double(abs(R(1))>0);
else
    p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
end

if p < ncolX
    % some columns need to be removed
    perm = perm(1:p); X = X(:,sort(perm));
end
end

function covariates = loadcovar(covarnames,var)
covariates = []; 
for i = 1:length(covarnames)
    if ismember(covarnames{i}, var.Properties.VariableNames) % is the variable name existing in the table?
        covariates = [covariates,var.(covarnames{i})];
        if any(isnan(covariates(:,end))) % missing observations?
            covariates(:,end) = [];
            disp([covarnames{i},' missing covariate removed'])
        end
    end
end
end