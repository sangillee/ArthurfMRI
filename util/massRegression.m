% massRegression
% function to run massive number of regression with same X but with different Y
% 2021-02-06: created
% 2021-02-14: renamed to massRegression, added additional measure calculations
% 2021-02-16: incorporated mean centering
% 2021-10-04: including the residual output
% 2022-03-22: incorporated safe calculation of zstat by always using negative tstat to calculate p-values

function [coef,outcome] = massRegression(X,Y,measure)
% if any(isnan(X(:))) % for debugging
%     keyboard
% end
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

switch measure
    case 'nothing'
        outcome = [];
    case 'R2'
        yhat = X*coef;
        outcome = diag(corr(yhat,Y))'.^2;
    case 'fstat'
        yhat = X*coef;
        R2 = diag(corr(yhat,Y))'.^2;
        outcome = (R2.*(n-p))./((1-R2).*(p-1));
    case 'zfstat'
        yhat = X*coef;
        R2 = diag(corr(yhat,Y))'.^2;
        F = (R2.*(n-p))./((1-R2).*(p-1));
        p = fcdf(F,p-1,n-p,'upper');
        outcome = -norminv(p);
    case 'tstat'
        RI = R\eye(p);
        rmse = vecnorm(Y-X*coef,2,1)./sqrt(max(0,n-p));
        se = zeros(ncolX,size(Y,2));
        se(perm,:) = repmat(rmse,p,1).*repmat(sqrt(sum(abs(RI).^2,2)),1,size(Y,2));
        outcome = coef./se;
    case 'zstat'
        RI = R\eye(p);
        rmse = vecnorm(Y-X*coef,2,1)./sqrt(max(0,n-p));
        se = zeros(ncolX,size(Y,2));
        se(perm,:) = repmat(rmse,p,1).*repmat(sqrt(sum(abs(RI).^2,2)),1,size(Y,2));
        tstat = coef./se;
        poststat = tstat>0; tstat(poststat) = -tstat(poststat); % calculating probabilities close to 0 is numerically safer than those close to 1
        pval = tcdf(tstat,n-p);
        outcome = norminv(pval);
        outcome(poststat) = -outcome(poststat);
    case 'resid'
        outcome = Y-X*coef;
    otherwise
        error('unknown output statistic requested')
end
end