function create_aal_mask(size,name)
aalimg = load_nii('aal.nii.gz');  % this is 181 217 181
aalmask = double(aalimg.img ~= 0 & aalimg.img < 91); % cerebrum only
aalmask = imresize3(aalmask,size); % fmriprep dimension is 97 115 97
aalmask = 1.*(aalmask>0.5);
Niftisave(name,aalmask)
end