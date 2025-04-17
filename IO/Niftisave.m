% function for easily saving matrix data into a nifti file
% smart-guessing of image resolution depending on matrix size

function info = Niftisave(data,name)
    data(isnan(data)) = 0; % nans are removed

    % first, what is the current dimension of the image?
    dims = size(data);
    switch length(dims)
        case 2 % this is either a 1D or 2D image
            if dims(1) == 1 % this is a row vector
                if dims(2) == 1; error('no scalars please'); end
                data = data'; dims = size(data); nT = 1; % convert to column vector
            elseif dims(2) == 1 % this is a column vector
                nT = 1;
            else % this is a matrix
                if dims(2) > dims(1) % assume that number of observations are less than number of voxels and convert all to column vectors
                    data = data'; dims = size(data);
                end
                nT = dims(2);
            end
            data = sizebank1D(data,dims(1),nT); dims = size(data); % convert to 3D/4D image
        case 3 % this is a 3d image
            nT = 1;
        case 4 % this is a 4d image
            nT = dims(4);
    end

    info = getTemplate(dims(1:3));
    
    info.Datatype = class(data);
    if nT>1 % adjust header if necessary
        info.raw.dim(1) = 4;
        info.ImageSize(4) = nT;
        info.raw.dim(5) = nT;
        info.raw.pixdim(5) = 2; % assuming 2 second TR
        info.PixelDimensions(4) = 2;
    end

    if nargin == 2
        niftiwrite(data,name,info,'Compressed',true)
    end
end

% a bank for querying the size of a vectorized brain image and loading appropriate images
function outimg = sizebank1D(imagedata,nvox,nT)
    switch nvox
        case 1082035 % standard fmripep 2mm dimension image, but vectorized
            mask = ones(97,115,97);
        case 902629 % standard fsl 2mm dimension image, but vectorized
            mask = ones(91,109,91);
        case 271633 % standard fsl 3mm dimension image, but vectorized
            mask = ones(61,73,61);
        case 235840 % vectorized fmriprep 2mm brain mask
            load('mask_brain_fmriprep_2mm.mat','mask')
        case 31538 % vectorized fmriprep 4mm brain mask
            load('mask_brain_fmriprep_4mm.mat','mask')
        case 231360 % vectorized fsl 2mm brain mask
            load('mask_brain_fsl_2mm.mat','mask')
        case 80321 % vectorized fsl 3mm brain mask
            load('mask_fsl_3mm.mat','mask')
        case 170074 % vectorized fmriprep 2mm Cerebrum Mask (smoothed by dil/ero)
            load('mask_aalcerebrum_fmriprep_2mm.mat','mask')
        case 24301 % vectorized fmriprep 4mm Cerebrum Mask (re-sampled from 2mm)
            load('mask_aalcerebrum_fmriprep_4mm.mat','mask')
        case 176230 % vectorized fsl 2mm Cerebrum Mask (smoothed by dil/ero)
            load('mask_aalcerebrum_fsl_2mm.mat','mask')
        case 52288 % vectorized fsl 3mm Cerebrum Mask (smoothed by dil/ero)
            load('mask_aalcerebrum_fsl_3mm.mat','mask')
        otherwise
            error('unexpected number of voxels')
    end
    outimg = zeros(length(mask(:)),nT);
    outimg(mask(:)==1,:) = imagedata;
    outimg = reshape(outimg,[size(mask),nT]);
end

function info = getTemplate(sizeinfo)
    if isequaln(sizeinfo,[97,115,97]) % fmriprep 2mm image
        load('template_fmriprep_2mm.mat','info')
    elseif isequaln(sizeinfo,[49,58,49]) % fmriprep 4mm image
        load('template_fmriprep_4mm.mat','info')
    elseif isequaln(sizeinfo,[91,109,91]) % FSL 2mm image
        load('template_fsl_2mm.mat','info')
    elseif isequaln(sizeinfo,[61,73,61]) % FSL 3mm image
        load('template_fsl_3mm.mat','info')
    else
        error('unknown template')
    end
end