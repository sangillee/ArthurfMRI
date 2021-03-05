% function for porting image in nativemask space into another space (usually mnispace) using a new mask

function mniimg = portimg(nativebeta,nativemask,mnimask)
mniimg = zeros(size(nativebeta,1),sum(mnimask(:)));
mnioverlap = nativemask(mnimask(:)==1);
nativeoverlap = mnimask(nativemask(:)==1);
mniimg(:,mnioverlap==1) = nativebeta(:,nativeoverlap==1);
end