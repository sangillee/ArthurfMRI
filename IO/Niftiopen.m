function img = Niftiopen(name)
    info = niftiinfo(name);
    img = niftiread(name);
    img = info.raw.scl_inter + info.raw.scl_slope*double(img); % apply scaling as prescribed in the header
end