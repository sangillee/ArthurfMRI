function [clusterinvp,rawinvpmap,statmap] = myPermute(brain,mask,initp,covar,covar0)
% Function for performing various permutation testing on brain maps
% inputs 
%   brain:  n-by-v matrix of brain images
%   mask:   3D brain mask for the input 'brain'
%   initp:  initial thresholding p-value. Default is 0.001;
%   covar:  (optional) a n-by-1 vector. If supplied, the function performs permutation testing for significant correlation with Y
%   covar0: (optional) a scalar. If supplied, the function calculates the predicted brain map when Y = Y0 using regression. This predicted
%           map is permutation tested

% input handling
if nargin < 5; covar0 = []; end
if nargin < 4; covar = []; end
if nargin < 3; initp = 0.001; end
init_thresh = 1-initp;

% options
nrep = 5000;

% analysistype
if isempty(covar) % flip-sign t-tests
    myfunction = @(braindata,covariate,covariate0,shuffle) ttestinvp(braindata, mask,shuffle);
elseif isempty(covar0) % correlation test
    myfunction = @(braindata,covariate,covariate0,shuffle) corrinvp(braindata,covariate, mask,shuffle);
else % regression prediction test
    myfunction = @(braindata,covariate,covariate0,shuffle) regpredinvp(braindata,covariate,covariate0, mask,shuffle);
end

% get original test p-value and clusters formed
[rawinvpmap,statmap] = myfunction(brain,covar,covar0,0);
origCC = bwconncomp(rawinvpmap>init_thresh);

% space allocation
maxstats = nan(nrep,1);
clusterinvp = zeros(size(mask));

if origCC.NumObjects > 0 % there's at least one cluster at initial threshold
    % do permutation with cluster-size thresholding
    parfor i = 1:nrep
        disp(i)
        [sampleinvpmap,samplestatmap] = myfunction(brain,covar,covar0,1); % resulting correlation p-values from shuffling
        CC = bwconncomp(sampleinvpmap>init_thresh); % identify clusters based on initial thresholding
        if isempty(CC.PixelIdxList)
            maxstats(i) = 0;
        else
            statstat = nan(length(CC.PixelIdxList),1);
            for j = 1:length(CC.PixelIdxList)
                statstat(j) = sum(samplestatmap(CC.PixelIdxList{j}));
            end
            maxstats(i) = max(statstat);
            % maxstats(i) = max(cellfun(@length,CC.PixelIdxList)); cluster size thresholding
        end
    end
    
    % assigning inverse p-values to clusters
    for i = 1:origCC.NumObjects
        clusterind = origCC.PixelIdxList{i};
        % clusterstat = length(clusterind); % cluster size thresholding
        clusterstat = sum(statmap(clusterind));
        invp = sum(maxstats<=clusterstat)./nrep;
        clusterinvp(clusterind) = invp;
    end
end
end

function [invpmap,rmap] = corrinvp(brain,covar,mask,shuffle)
if shuffle == 1
    n = length(covar);
    shuffling = randsample(n,n); % random shuffling of pairs
    [r,p] = corr(brain,covar(shuffling));
else
    [r,p] = corr(brain,covar);
end
invpmap = mask;
invpmap(mask(:)==1) = 1-p; % inverse p-value
rmap = mask;
rmap(mask(:)==1) = r;
end

function [invpmap,tmap] = ttestinvp(brain,mask,shuffle)
if shuffle == 1
    shuffle = 2.*(randi(2,size(brain,1),1)-1.5);
    [~,p,~,stats] = ttest(brain.*shuffle);
else
    [~,p,~,stats] = ttest(brain);
end
invpmap = mask;
invpmap(mask(:)==1) = 1-p; % inverse p-value
tmap = mask;
tmap(mask(:)==1) = stats.tstat;
end