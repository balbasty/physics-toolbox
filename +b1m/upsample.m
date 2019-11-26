function s = upsample(s, lat, bnd)
% Upsample the sensitivity fields to a given lattice.
%
% FORMAT s = b1m.upsample(s, lat, [bnd])
%
% s   - Input volume
% lat - Output lattice dimensions
% bnd - Boundary condition about each dimension (0=circulant/[1]=m)
%
% This upsampling is obtained by performing a discrete transform (DCT/FFT
% for Mirror/Circulant conditions), zero-padding, and
% performing the inverse transform.
%
% This function relies on Matlab's signal processing toolbox, and might
% need a recent version of Matlab (>= 2016).
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Initial and final grids
lat0 = size(s);
ndim = max(numel(lat0), numel(lat));
lat0 = utils.pad(lat0, [0 ndim-numel(lat0)], 1, 'post');
lat  = utils.pad(lat,  [0 ndim-numel(lat)],  Inf, 'post');
lat(~isfinite(lat)) = lat0(~isfinite(lat));
if nargin < 3
    bnd = 1;
end
bnd  = utils.pad(bnd,  [0 ndim-numel(bnd)],  'replicate', 'post');

for i=1:numel(lat)
    if isfinite(lat(i)) && (lat(i) ~= lat0(i))
        switch bnd(i)
            case 0 % Circulant
                s      = utils.fft(s, i);
                pad    = zeros(1,ndim);
                pad(i) = floor((lat(i)-lat0(i))/2);
                s      = utils.pad(s, pad, 0, 'both');
                if mod(lat(i)-lat0(i),2)
                    % odd padding
                    pad(i) = 1;
                    if mod(lat0(i),2), dir = 'pre';
                    else,              dir = 'post'; end
                    s = utils.pad(s, pad, 0, dir);
                end
                s = utils.ifft(s, i);
            case 1 % Mirror
                s = sqrt(lat(i)/lat0(i)) * idct(dct(s,[],i),lat(i),i);
        end
    end
end