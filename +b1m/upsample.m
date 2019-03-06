function s = upsample(s, lat)
% Upsample the sensitivity fields to a given lattice.
%
% FORMAT s = b1m.upsample(s, lat)
%
% This upsampling is obtained by performing a discrete cosine transform (in
% order to keep Neumann's boundary conditions), zero-padding, and
% performing the inverse DCT.
%
% This function relies on Matlab's signal processing toolbox, and might
% need a recent version of Matlab (>= 2016).
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

lat0 = size(s);
lat0 = padarray(lat0, [0 max(0,numel(lat)-numel(lat0))], 1, 'post');

for i=1:numel(lat)
    if isfinite(lat(i)) && (lat(i) ~= lat0(i))
        s = sqrt(lat(i)/lat0(i)) * idct(dct(s,[],i),lat(i),i);
    end
end