% searchlight_getNeighbors
% Arthur Lee
% Written 1/31/2021
% Revised 7/28/2022
% function that gets the neighboring voxels of each voxel inside a 3d mask image
% Usage: neighborMatrix = getSearchlightNeighbors(mask3d,dist)
% Output: output matrix that carries, for each voxel inside mask3d, a row of voxel indices that belong inside that voxel's searchlight
% Input: mask3d (3d binary mask image), dist (radius of searchlight sphere in voxels)

function neighborList = searchlight_getNeighbors(mask3d,dist)
v = sum(mask3d(:)); [xsize,ysize,zsize] = size(mask3d);
offsets = getOffsets('dist',dist); % create searchlight voxel coordinate offset list
neighborList = cell(v,1); % output matrix
[x,y,z] = ind2sub([xsize,ysize,zsize],find(mask3d(:)==1)); coord = [x,y,z]; % x y z coordinates of mask voxels
betaind = cumsum(mask3d(:)); betaind(mask3d(:)==0) = nan; % indice from 1 to v for all voxels in the mask
for i = 1:v
    xyzind = [coord(i,:);coord(i,:) + offsets]; %xyz coordinates of potential searchlight voxels. center voxel plus neighbors
    
    % find neighbors that are outside the image bounds
    outsidebounds = any(xyzind<1,2) | xyzind(:,1)>xsize | xyzind(:,2)>ysize | xyzind(:,3)>zsize;
    xyzind(outsidebounds,:) = [];
    
    % convert back to linearInd
    searchlightind = sub2ind([xsize,ysize,zsize], xyzind(:,1),xyzind(:,2),xyzind(:,3));
    neighborList{i} = remnan(betaind(searchlightind)); % remove neighbors that are not inside the mask3d
end
end

function nonnan = remnan(input)
nonnan = input(~isnan(input));
end