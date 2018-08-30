% -------------------------------------------------------------------------
% Parameters
vs    = [1 1 1];    % > Voxel size
itmax = 1000;       % > Maximum number of optim iterations
prm   = [0 0 1E7];  % > Sensitivity fields regularisation 
                    %   ['abs value' 'membrane E' 'bending E']

% Note: The bending energy is expressed in mm^(-4)
%       A way of setting this paramerer can thus be [Expected FWHM in mm]^4
%           FWHM = 50mm => prm = 6.25*1E6

% -------------------------------------------------------------------------
% Phantom data
path = fileparts(which('test_multicoil_brain'));
load(fullfile(path,'phantom','phantom_data.mat'));
kData = reshape(kData, [size(kData,1) size(kData,2) 1 size(kData,3)]);

% -------------------------------------------------------------------------
% Downsample
downsample = true;
dim        = [32 32 1];
if downsample
    dim0   = [size(kData,1) size(kData,2) size(kData,3)];
    factor = dim0./dim;
    center = ceil(dim0/2);
    imin   = floor(center-dim/2)+1;
    imax   = floor(center+dim/2);
    kData  = kData(imin(1):imax(1),imin(2):imax(2),imin(3):imax(3),:);
    vs = vs .* factor;
end

% -------------------------------------------------------------------------
% Fourier transform
x = double(kData); % clear kData
x = fftshift(x, 1);
x = fftshift(x, 2);
x = ifft(x, [], 1);
x = ifft(x, [], 2);
x = fftshift(x, 1);
x = fftshift(x, 2);
x = cat(5, real(x), imag(x));
x = single(x);

% -------------------------------------------------------------------------
% Select number of channels
N = size(x,4);      % < e.g., 2
x = x(:,:,:,1:N,:);

% -------------------------------------------------------------------------
% Allocate unknown arrays
s   = zeros(size(x));                                   % < sensitivity
rho = zeros(size(x,1),size(x,2),size(x,3),1,size(x,5)); % < mean

% -------------------------------------------------------------------------
% Noise estimates
% > We start with some sort of initial estimate based on a very
%   simple independant Gaussian (mixture) fitting
[C,A] = multicoil_init_cov(x); % < C = covariance / A = inv(C) = precision

% -------------------------------------------------------------------------
% Main
[s,rho,A,ll] = multicoil_infer(x, s, rho, A, prm, vs, itmax, 2);