% massRegression
% function to run massive number of regression with same X but with different Y
% 2021-02-06: created
% 2021-02-14: renamed to massRegression, added additional measure calculations
% 2021-02-16: incorporated mean centering
% 2021-10-04: including the residual output
% 2022-03-22: incorporated safe calculation of zstat by always using negative tstat to calculate p-values

function [coef,stats] = massRegression(X,Y,measure)
assert(~any(isnan(X(:))),'nan found in X matrix')
assert(~any(isnan(Y(:))),'nan found in Y matrix')

% mean centering to remove the need for intercept
X = X-mean(X); Y = Y-mean(Y);

% Use the rank-revealing QR to remove dependent columns of X.
[n,ncolX] = size(X); [Q,R,perm] = qr(X,0);
if isempty(R)
    p = 0;
elseif isvector(R)
    p = double(abs(R(1))>0);
else
    p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
end
if p < ncolX
    warning(message('stats:regress:RankDefDesignMat'));
    R = R(1:p,1:p); Q = Q(:,1:p); perm = perm(1:p);
end
coef = zeros(ncolX,size(Y,2));coef(perm,:) = R \ (Q'*Y);

if strcmp(measure,'nothing')
    stats = [];
else
    yhat = X*coef;
    switch measure
        case 'R2'
            stats = (sum(zscore(yhat).*zscore(Y))./ (n-1)).^2; %diag(corr(yhat,Y))'.^2;
        case 'fstat'
            R2 = diag(corr(yhat,Y))'.^2;
            stats = (R2.*(n-p))./((1-R2).*(p-1));
        case 'zfstat'
            R2 = diag(corr(yhat,Y))'.^2;
            F = (R2.*(n-p))./((1-R2).*(p-1));
            pval = fcdf(F,p-1,n-p,'upper');
            stats = -norminv(pval);
        case 'tstat'
            stats = getT(R,p,Y,yhat,n,ncolX,perm,coef);
        case 'zstat'
            tstat = getT(R,p,Y,yhat,n,ncolX,perm,coef);
            poststat = tstat>0; tstat(poststat) = -tstat(poststat); % calculating probabilities close to 0 is numerically safer than those close to 1
            pval = tcdf(tstat,n-p-1);
            stats = norminv(pval);
            stats(poststat) = -stats(poststat);
        case 'resid'
            stats = Y-yhat;
        otherwise
            error('unknown output statistic requested')
    end
end
end

function tstat = getT(R,p,Y,yhat,n,ncolX,perm,coef)
rmse = vecnorm(Y-yhat,2,1)./sqrt(max(0,n-p));
se = zeros(ncolX,size(Y,2));
se(perm,:) = repmat(rmse,p,1).*repmat(sqrt(sum(abs(R\eye(p)).^2,2)),1,size(Y,2));
tstat = coef./se;
end