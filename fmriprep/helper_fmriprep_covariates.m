% helper_fmriprep_covariates
% function to take care of loading images, masks, and covariates
% 2021-02-06: created
% 2021-02-16: renamed to helper_fmriprep_Covariates from helper_fmriprep_defaultCovariates
%             now accepts option structure to identify which covariates to include

function covariates = helper_fmriprep_covariates(var,opts)
if nargin <2
    opts = helper_fmriprep_regoptions; % default options
    opts = opts.covar;
end
errorchecking(opts) % option struct error checking
covariates = []; % matrix for covariates

% cosine components
covarnames = cell(1,0);
for i = 1:opts.cosines
    covarnames{end+1} = ['cosine',sprintf('%02d',i-1)]; % cosine components
end
covariates = [covariates,loadcovar(covarnames,var)];

% a_comp_cor components
covarnames = cell(1,0);
for i = 1:opts.a_comp_cor
    covarnames{end+1} = ['a_comp_cor_',sprintf('%02d',i-1)]; % adding top CSF+WM components.
end
covariates = [covariates,loadcovar(covarnames,var)];

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
    case 6
        motion = loadcovar({'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'},var);
    case 12
        motion = loadcovar({'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z','trans_x_power2','trans_y_power2','trans_z_power2','rot_x_power2','rot_y_power2','rot_z_power2'},var);
    case 24
        motion = loadcovar({'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z','trans_x_power2','trans_y_power2','trans_z_power2','rot_x_power2','rot_y_power2','rot_z_power2'},var);
        motion = [motion,[mean(motion);motion(1:end-1,:)]]; % t-1 motion
        motion(:,end+1) = [1;zeros(size(motion,1)-1,1)]; % first volume index
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

function errorchecking(opts)
assert(isstruct(opts),'covariate option structure is not a struct')
expectation = {'cosines';'a_comp_cor';'csf';'white_matter';'motionparams'};
names = fieldnames(opts);
assert(isequaln(names,expectation),'unexpected fields in covariate option')

% all options must be numeric, non-negative integers
for i = 1:length(names)
    assert(isnumeric(opts.(names{i})),'all covariate option values must be numeric')
    assert(opts.(names{i})>=0,'all covariate option values must be non-negative')
    assert(floor(opts.(names{i}))==opts.(names{i}),'all covariate option values must be integers')
end

% additional checks
assert(ismember(opts.csf,[0,1]),'csf option should be either 0 or 1')
assert(ismember(opts.white_matter,[0,1]),'white_matter option should be either 0 or 1')
assert(ismember(opts.motionparams,[0,6,12,24]),'motion param options should be one of 0, 6, 12, or 24')
end