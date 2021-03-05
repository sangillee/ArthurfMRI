% created Feb 26th 2021
% function to subsample a 3D (or 4D) image into half resolution

function betas = subsamp2(coef,origmask,downmask)
origmask = double(origmask);
betas = nan(size(coef,1),sum(downmask(:)));
for i = 1:size(coef,1)
    img = origmask;
    img(origmask(:)==1) = coef(i,:); % put beta into 3d image
    img = imresize3(img,0.5,'method','linear'); % linear down sampling interpolation
    betas(i,:) = img(downmask(:)==1);
end
end