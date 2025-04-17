function myPermute_auto(filename,brain,mask,initp,covar,stattype,correctiontype,corrtype)
if nargin < 8; corrtype = 'Spearman'; end   % spearman
if nargin < 7; correctiontype = 'FDR'; end  % FDR
if nargin < 6; stattype = 'mass'; end       % cluster mass thresholding
if nargin < 5; covar = []; end
if nargin < 4; initp = 0.001; end

[clusterinvp1,clusterinvp2,statmap] = myPermute(brain,double(mask),initp,covar,stattype,correctiontype,corrtype);
maxinvp = max(clusterinvp1,clusterinvp2);

initpstring = num2str(initp);
map = cat(4,(maxinvp>0.95).*statmap,statmap,maxinvp);
Niftisave(map,[filename,'_',stattype,'_',correctiontype,'_',initpstring(3:end)]);
end