function image = smoothfMRI(image,FWHM,voxsize,mask)
if nargin < 4
    mask = [];
end
if nargin < 3
    voxsize = 2;
end

if isempty(mask)
    for i = 1:size(image,4)
        image(:,:,:,i) = imgaussfilt3(image(:,:,:,i),FWHM/2.355/voxsize);
    end
else
    for i = 1:size(image,1)
        tempimg = double(mask);
        tempimg(mask(:)==1) = image(i,:);
        tempimg = imgaussfilt3(tempimg,FWHM/2.355/voxsize);
        image(i,:) = tempimg(mask(:)==1);
    end
end
end