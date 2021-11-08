function reg = simBOLD(TR,nvol,fslEV,type)
if nargin<4
    type = 'spm';
end
dt = 0.01;
boxcar_y = BOXCAR(TR*nvol,fslEV(:,1),fslEV(:,1)+fslEV(:,2),fslEV(:,3),dt);
if strcmp(type,'fsl')
    signal = conv(local_gamma_hrf(dt),boxcar_y);
    cutind = (length(signal)-length(boxcar_y));
    signal = signal(1:end-cutind);
elseif strcmp(type,'spm')
    signal = conv(local_spm_hrf(dt),boxcar_y);
    cutind = (length(signal)-length(boxcar_y));
    signal = signal(1:end-cutind);
elseif strcmp(type,'box')
    signal = boxcar_y;
else
    error('unknown HRF requested')
end

reg = signal(1:(TR/dt):(end-1));
% if you want to see the results:
if 0
    subplot(3,1,1)
    plot(boxcar_y)
    subplot(3,1,2)
    plot(signal)
    subplot(3,1,3)
    plot(reg)
end
end

function y = BOXCAR(dur,onset,offset,mod,dt)
y = repmat(0:dt:dur,length(onset),1); ind = y<onset | y>offset;
y(ind) = 0; y(~ind) = 1;
mod(isnan(mod)) = 0;
y = y'*mod;
end

function hrf = local_spm_hrf(RT)
% p    - parameters of the response function (two Gamma functions)
%                                                           defaults
%                                                          {seconds}
%        p(1) - delay of response (relative to onset)          6
%        p(2) - delay of undershoot (relative to onset)       16
%        p(3) - dispersion of response                         1
%        p(4) - dispersion of undershoot                       1
%        p(5) - ratio of response to undershoot                6
%        p(6) - onset {seconds}                                0
%        p(7) - length of kernel {seconds}                    32

p = [6 16 1 1 6 0 32];
fMRI_T = 16; % microtime resolution
dt  = RT/fMRI_T;
u   = (0:ceil(p(7)/dt)) - p(6)/dt;
hrf = spm_Gpdf(u,p(1)/p(3),dt/p(3)) - spm_Gpdf(u,p(2)/p(4),dt/p(4))/p(5);
hrf = hrf((0:floor(p(7)/RT))*fMRI_T + 1);
hrf = hrf'/sum(hrf);
end

function hrf = local_gamma_hrf(RT)
t = 0:RT:32;
meanlag = 6; stddev = 3;
a = (meanlag/stddev)^2;
b = meanlag/(stddev^2);
hrf = pdf_gamma(t,a,b) * 0.0063;
end

function pdfx = pdf_gamma(x,a,b)
nx = length(x);
pdfx = zeros(nx,1);
ind = find(x>0);
pdfx(ind) = (b.^2) .* (x(ind).^(a-1)) .* exp(-b.*x(ind)) ./ gamma(a);
end