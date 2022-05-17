% [outcome,coef,X,covariates] = fitfMRI(REG,Y,covariates,TR,HP,measure)
% function for regressing Y against the convolved version of REG and covariates
% REG is a cell array where each cell contains a three-column FSL-style matrix (first column: onset, second column: duration, third column: modulation)
% Each cell will create a expected HRF signal and be entered as regressors
% covariates are variables that you want to use as regressors without convolving
% HP is a switch for FSL-style High pass filtering on Y and X (but not covariates). HP is the cutoff seconds (i.e., 50 sec), or 0 if turned off.
% Y may be a matrix or a vector. If it's a matrix, each column of Y will be regressed using the same predictor variables
% measure is the statistic that you want. can be MSE, R2, adjR2
% calls on simBOLD

% 2022-03-29: Added option for temporal derivative, which can be used by having a fourth column on the FSL-style matrix for REG

function [coef,stats] = fitfMRI(REG,Y,covariates,TR,HP,measure,hrf,drop)
if nargin < 8; drop = 0; end
if nargin < 7; hrf = 'spm'; end
if nargin < 6; measure = 'nothing'; end
if nargin < 5; HP = 0; end
nvol = size(Y,1);

% create regressors
if ~isempty(REG)
    X = nan(nvol,length(REG));
    tempderiv = nan(nvol,length(REG));
    tempcounter = 0;
    for i = 1:size(X,2)
        if size(REG{i},2) == 4 % temporal derivative requested
            tempcounter = tempcounter+1;
            [X(:,i),tempderiv(:,tempcounter)] = simBOLD(TR,nvol,REG{i}(:,1:3),hrf);
        else
            X(:,i) = simBOLD(TR,nvol,REG{i}(:,1:3),hrf);
        end
    end
    X = [X,tempderiv(:,1:tempcounter),covariates];
else
    X = covariates;
end

if drop > 0
    X(1:drop,:) = [];
    Y(1:drop,:) = [];
end

% if high-pass filtering is used
if HP > 0
    F = HPfilter(nvol-drop,HP,TR); % filter calculation
    Y = F*Y; X = F*X;
end

% run regression
[coef,stats] = massRegression(X,Y,measure);
end

function F = HPfilter(nvol,cutoffseconds,TR)
sigN2=((cutoffseconds/TR)/sqrt(2))^2;
K=toeplitz(1/sqrt(2*pi*sigN2)*exp(-(0:(nvol-1)).^2/(2*sigN2)));
coder.extrinsic('spdiags');
K=spdiags(1./sum(K,2), 0, nvol,nvol)*K;
H = zeros(nvol,nvol); % Smoothing matrix, s.t. H*y is smooth line
Q = [ones(nvol,1) (1:nvol)'];
for k = 1:nvol
    W = diag(K(k,:));
    Hat = Q*pinv(W*Q)*W;
    H(k,:) = Hat(k,:);
end
F=eye(nvol)-H;
end