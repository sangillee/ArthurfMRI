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