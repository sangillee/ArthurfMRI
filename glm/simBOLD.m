% A function for simulating BOLD response based on a fsl-style three-column regressor (first column: onset (sec), second column: duration (sec), third column: modulation value)
%
% inputs:
%   TR: repetition time, in seconds
%   nvol: total number of volumes (time points) requested
%   fslEV: a n-by-3 column matrix in fsl custom regressor format
%   type: convolution kernel to use. 'fsl' for single gamma, 'spm' for double gamma, and 'box' for boxcar
%
% outputs:
%   reg: generated BOLD time series
%   tD: temporal derivative of the generated reg. only calculated when requested
function [reg,tD] = simBOLD(TR,nvol,fslEV,type)
    if nargin<4
        type = 'fsl'; % default HRF type of single gamma
    end
    dt = 0.05; % time resolution

    % generate high-temporal-resolution regressors
    boxcar_y = BOXCAR(TR*nvol,fslEV(:,1),fslEV(:,1)+fslEV(:,2),fslEV(:,3),dt);
    if strcmp(type,'fsl')
        signal = conv(local_fsl_hrf(dt),boxcar_y);
    elseif strcmp(type,'spm')
        signal = conv(local_spm_hrf(dt),boxcar_y);
    elseif strcmp(type,'box')
        signal = boxcar_y;
    else
        error('unknown HRF requested')
    end

    % down-sampling to match middle of each TR
    reg = signal(round(((1:nvol)-0.5)*TR/dt)+1);

    % temporal derivative calculation, if requested
    if nargout > 1
        tD = [reg(2)-reg(1);0.5.*(reg(3:end)-reg(1:end-2));reg(end)-reg(end-1)];
        temp = [reg,ones(size(reg))];
        tD = tD - temp*mldivide(temp,tD); % orthogonalize
    end
end

function y = BOXCAR(dur,onset,offset,mod,dt)
    y = repmat(0:dt:dur,length(onset),1); ind = y<onset | y>=offset;
    y(ind) = 0; y(~ind) = 1;
    mod(isnan(mod)) = 0;
    y = y'*mod;
end

function hrf = local_fsl_hrf(DT)
    t = 0:DT:(21-DT);
    meanlag = 6; stddev = 3;
    a = (meanlag/stddev)^2;
    b = (stddev^2)/meanlag;
    hrf = gampdf(t,a,b)*DT;%002405;
    hrf = hrf./sum(hrf);
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
    hrf = localspm_Gpdf(u,p(1)/p(3),dt/p(3)) - localspm_Gpdf(u,p(2)/p(4),dt/p(4))/p(5);
    hrf = hrf((0:floor(p(7)/RT))*fMRI_T + 1);
    hrf = hrf'/sum(hrf);
end

function f = localspm_Gpdf(x,h,l)
    if nargin<3, error('Insufficient arguments'), end

    ad = [ndims(x);ndims(h);ndims(l)];
    rd = max(ad);
    as = [[size(x),ones(1,rd-ad(1))];...
        [size(h),ones(1,rd-ad(2))];...
        [size(l),ones(1,rd-ad(3))]];
    rs = max(as);
    xa = prod(as,2)>1;
    if sum(xa)>1 && any(any(diff(as(xa,:)),1))
        error('non-scalar args must match in size');
    end

    f = zeros(rs);

    md = ( ones(size(x))  &  h>0  &  l>0 );
    if any(~md(:))
        f(~md) = NaN;
        warning('Returning NaN for out of range arguments');
    end

    ml = ( md  &  x==0  &  h<1 );
    f(ml) = Inf;
    ml = ( md  &  x==0  &  h==1 ); if xa(3), mll=ml; else mll=1; end
    f(ml) = l(mll);

    %-Compute where defined and x>0
    Q  = find( md  &  x>0 );
    if isempty(Q), return, end
    if xa(1), Qx=Q; else Qx=1; end
    if xa(2), Qh=Q; else Qh=1; end
    if xa(3), Ql=Q; else Ql=1; end

    %-Compute
    f(Q) = exp( (h(Qh)-1).*log(x(Qx)) +h(Qh).*log(l(Ql)) - l(Ql).*x(Qx)...
        -gammaln(h(Qh)) );
end