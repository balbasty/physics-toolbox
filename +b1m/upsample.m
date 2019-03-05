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

lat0 = [size(s) 1];
lat0 = lat0(1:3);

s = dct(dct(dct(s, [], 1), [], 2), [], 3);
s = idct(idct(idct(s, lat(1), 1), lat(2), 2), lat(3), 3);
s = s * sqrt(prod(lat)/prod(lat0));
