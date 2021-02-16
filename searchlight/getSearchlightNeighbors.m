% getSearchlightNeighbors
% Arthur Lee
% Written 1/31/2021
% function that gets the neighboring voxels of each voxel inside a 3d mask image
% Usage: neighborMatrix = getSearchlightNeighbors(mask3d,dist)
% Output: output matrix that carries, for each voxel inside mask3d, a column of voxel indices that belong inside that voxel's searchlight
% Input: mask3d (3d binary mask image), dist (radius of searchlight sphere in voxels)

function neighborMatrix = getSearchlightNeighbors(mask3d,dist)
offsets = createXYZoffsets(dist); % create searchlight voxel coordinate offset list
neighborMatrix = nan(size(offsets,1),sum(mask3d(:))); % output matrix
linearind = find(mask3d(:)==1); % linear indice of alive voxels in the 3d image
betaind = cumsum(mask3d(:)); betaind(mask3d(:)==0) = nan; % indice from 1 to v for all voxels in the mask
for i = 1:length(linearind)
    i
    [x,y,z] = ind2sub(size(mask3d),linearind(i)); % convert linear indice to x y z coordinates
    xyzind = [x+offsets(:,1),y+offsets(:,2),z+offsets(:,3)]; %xyz coordinates of potential searchlight voxels
    % find neighbors that are outside the image bounds
    xyzind(any(xyzind<1,2),:) = [];
    xyzind(xyzind(:,1)>size(mask3d,1),:) = [];
    xyzind(xyzind(:,2)>size(mask3d,2),:) = [];
    xyzind(xyzind(:,3)>size(mask3d,3),:) = [];
    % convert back to linearInd
    searchlightind = sub2ind(size(mask3d), xyzind(:,1),xyzind(:,2),xyzind(:,3));
    % remove neighbors that are not inside the mask3d
    searchlightind = remnan(betaind(searchlightind));
    neighborMatrix(1:length(searchlightind),i) = searchlightind;
end
end

function nonnan = remnan(input)
nonnan = input(~isnan(input));
end

function offsets = createXYZoffsets(dist)
% create a triplet matrix of coordinate offsets from voxel of interest within a given distance
rounddist = ceil(dist); % rounding up the distance measure to nearest whole number
[xoffset,yoffset,zoffset] = ndgrid(-rounddist:rounddist,-rounddist:rounddist,-rounddist:rounddist); % candidate offsets near the voxel of interest
offsets = [xoffset(:),yoffset(:),zoffset(:)]; % triplet form of offsets
voxdist = sqrt(sum(offsets.^2,2)); % euclidean distance
offsets(voxdist>dist,:) = []; % remove far ones
end