% helper_fmriprep_regoptions
% function to make option struct to use in regression
% 2021-02-16: created with following default options:
% % opts = struct;
% % opts.FWHM = 0; % default no smoothing
% % opts.measure = 'nothing'; % default no statistic from regression
% % opts.HP = 0; % deafult no high-pass filtering
% % opts.HRF = 'spm'; % default hrf
% % opts.covar.cosines = 'all'; % cosine component covariates
% % opts.covar.a_comp_cor = 10; % a_comp_cor covariates
% % opts.covar.csf = 1; % average csf activity covariate
% % opts.covar.white_matter = 1; % average white matter activity covariates
% % opts.covar.motionparams = 24; % 0 for none, 6 for default, 12 for squareds, 24 for t-1 expansion

function opts = helper_fmriprep_regoptions
opts = struct;
opts.FWHM = 0; % default no smoothing
opts.measure = 'nothing'; % default no statistic from regression
opts.HP = 0; % deafult no high-pass filtering
opts.HRF = 'spm'; % default hrf
opts.covar.cosines = 9; % cosine component covariates
opts.covar.a_comp_cor = 10; % a_comp_cor covariates
opts.covar.csf = 1; % average csf activity covariate
opts.covar.white_matter = 1; % average white matter activity covariates
opts.covar.motionparams = 24; % 0 for none, 6 for default, 12 for squareds, 24 for t-1 expansion
end