% [outcome,coef,X,covariates] = fitfMRI(REG,Y,covariates,TR,HP,measure)
% function for regressing Y against the convolved version of REG and covariates
% REG is a cell array where each cell contains a three-column FSL-style matrix (first column: onset, second column: duration, third column: modulation)
% Each cell will create a expected HRF signal and be entered as regressors
% covariates are variables that you want to use as regressors without convolving
% HP is a switch for FSL-style High pass filtering on Y and X (but not covariates). HP is the cutoff seconds (i.e., 50 sec), or 0 if turned off.
% Y may be a matrix or a vector. If it's a matrix, each column of Y will be regressed using the same predictor variables
% measure is the statistic that you want. can be MSE, R2, adjR2
% calls on simBOLD

function [coef,outcome] = fitfMRI(REG,Y,covariates,TR,HP,measure,hrf)
if nargin < 7; hrf = 'spm'; end
if nargin < 6; measure = 'nothing'; end
nvol = size(Y,1);

% create regressors
X = nan(nvol,length(REG));
for i = 1:size(X,2)
    X(:,i) = simBOLD(TR,nvol,REG{i},hrf);
end

% if high-pass filtering is used
if HP > 0
    F = FSLHPfilter(nvol,HP,TR); % filter calculation
    Y = F*Y; X = F*X; covariates = F*covariates;
end

% run regression
[coef,outcome] = massRegression([X,covariates],Y,measure);
end