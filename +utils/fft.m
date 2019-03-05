function X = fft(X, DIMS)
% FORMAT Y = utils.fft(X, [DIMS])
% 
% Performs the discrete Fourier transform of the array X across the
% dimensions given in DIM.
% Conversely to Matlab's FFT, this function can work across multiple
% dimensions. Also, an IFFTSHIFT is performed prior to the FFT, and an
% FFTSHIFT is performed after.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

lat = size(X);
if nargin < 2
    DIMS = lat(find(lat>1, 1));
end
DIMS = DIMS(:)';

for d=DIMS
    X = fftshift(fft(ifftshift(X,d),[],d),d);
end