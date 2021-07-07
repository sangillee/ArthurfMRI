function R2out = varianceTiming(timingmatrix,Y,covariates,TR,HP,hrf)
if nargin < 6; hrf = 'spm'; end
if nargin < 5; HP = 0; end
nvol = size(Y,1);
numEpoch = size(timingmatrix,2);

R2out = nan(10*(numEpoch-1)+1,length(Y));
counter = 1;
for i = 1:.1:numEpoch
    % create regressors
    timing = genTime(timingmatrix,i);
    
    X = nan(nvol,size(timing,1));
    for j = 1:size(X,2)
        X(:,j) = simBOLD(TR,nvol,[timing(j), 0.1, 1],hrf);
    end
    X = [X,covariates];
    
    % if high-pass filtering is used
    if HP > 0
        F = HPfilter(nvol,HP,TR); % filter calculation
        Y = F*Y; X = F*X;
    end
    
    % run regression
    [~,R2out(counter,:)] = massRegression(X,Y,'R2'); counter = counter+1;
end
end

function timing = genTime(timingmatrix,interpolant)
startIND = floor(interpolant);
interpolant = interpolant - startIND;
startTime = timingmatrix(:,startIND);
if interpolant == 0
    timing = startTime;
else
    endTime = timingmatrix(:,startIND+1);
    timing = (1-interpolant).*startTime + interpolant.*endTime;
end

end