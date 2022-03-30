% helper_fmriprep_regoptions
% function to make option struct to use in regression
% 2021-02-16: created with default options
% 2022-03-27: created option checking function
%
% % default options
% % opts.FWHM = 0; % default no smoothing
% % opts.measure = 'nothing'; % no additional outputs by default
% % opts.HP = 0; % deafult no high-pass filtering
% % opts.HRF = 'spm'; % default hrf
% % opts.drop = 2; % default first few volumes to remove
% % opts.covar.cosines = 50; % maximum number of cosine component covariates. Surely, 50 would be enough even if the scan is long...
% % opts.covar.a_comp_cor = 10; % maximum number of a_comp_cor covariates
% % opts.covar.csf = 0; % average csf activity covariate
% % opts.covar.white_matter = 0; % average white matter activity covariates
% % opts.covar.motionparams = 24; % 0 for none, 6 for default, 12 for squareds, 24 for t-1 expansion

function opts = helper_fmriprep_regoptions(opts)
if nargin == 0
    % create new option struct
    opts = struct;
    opts.FWHM = 0;                  % default no smoothing
    opts.measure = 'nothing';       % no additional outputs by default
    opts.HP = 0;                    % deafult no high-pass filtering
    opts.HRF = 'spm';               % default hrf
    opts.drop = 2;                  % default first few volumes to remove
    opts.covar.cosines = 50;        % cosine component covariates
    opts.covar.a_comp_cor = 10;     % a_comp_cor covariates
    opts.covar.csf = 0;             % average csf activity covariate
    opts.covar.white_matter = 0;    % average white matter activity covariates
    opts.covar.motionparams = 24;   % 0 for none, 6 for default, 12 for squareds, 24 for t-1 expansion
else
    % perform error checking on the options
    assert(isstruct(opts),'option structure is not a struct')
    expectation = {'FWHM';'measure';'HP';'HRF';'drop';'covar'};
    names = fieldnames(opts);
    assert(isequaln(names,expectation),'unexpected fields in covariate option')
    
    % error checking for covariates
    assert(isstruct(opts.covar),'covariate option structure is not a struct')
    expectation = {'cosines';'a_comp_cor';'csf';'white_matter';'motionparams'};
    names = fieldnames(opts.covar);
    assert(isequaln(names,expectation),'unexpected fields in covariate option')
    
    % all options must be numeric, non-negative integers
    for i = 1:length(names)
        assert(isnumeric(opts.covar.(names{i})),'all covariate option values must be numeric')
        assert(opts.covar.(names{i})>=0,'all covariate option values must be non-negative')
        assert(floor(opts.covar.(names{i}))==opts.covar.(names{i}),'all covariate option values must be integers')
    end
    
    % additional checks
    assert(ismember(opts.covar.csf,[0,1]),'csf option should be either 0 or 1')
    assert(ismember(opts.covar.white_matter,[0,1]),'white_matter option should be either 0 or 1')
    assert(ismember(opts.covar.motionparams,[0,6,12,24,25]),'motion param options should be one of 0, 6, 12, or 24')
end
end