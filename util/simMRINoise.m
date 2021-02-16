function [noise,noiseinfo] = simMRINoise(nvol, sigma, n, AR)
% nvol = number of volumes
% sigma = noise standard deviation
% n = total number of noise vectors
% AR = autocorrelation coefficient
if nargin<4
    AR = 0.12;
end

corrmat = AR.^toeplitz(0:(nvol-1));
noiseinfo{1} = zeros(1,nvol);
noiseinfo{2} = (sigma^2).*corrmat;
noise = mvnrnd(noiseinfo{1},noiseinfo{2},n)';
end