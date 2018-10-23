% -------------------------------------------------------------------------
% Parameters
vs    = [8 8 8];    % > Voxel size
prm   = [1E7 1E7];  % > Sensitivity fields regularisation [LogMag Phase]

% Note: The bending energy is expressed in mm^(-4)
%       A way of setting this paramerer can thus be [Expected FWHM in mm]^4
%           FWHM = 50mm => prm = 6.25*1E6

% -------------------------------------------------------------------------
% Phantom data
path = fileparts(which('test_multicoil_brain'));
load(fullfile(path,'brain','complex_by_coil_images.mat'));
x = cat(5,real(rdata),imag(rdata));

% -------------------------------------------------------------------------
% Select number of channels
idx = randperm(size(x,4));
N = size(x,4);      % < e.g., 2
% N = 2;
x = x(:,:,:,idx(1:N),:);

% -------------------------------------------------------------------------
% Allocate unknown arrays

% -------------------------------------------------------------------------
% Noise estimates
% > We start with some sort of initial estimate based on a very
%   simple independant Gaussian (mixture) fitting
[C,A] = multicoil_init_cov(x); % < C = covariance / A = inv(C) = precision

% -------------------------------------------------------------------------
% Initialise on magnitude
if true
    mag = sqrt(sum(x.^2,5));
    s   = zeros(size(mag),'single');                     % < sensitivity
    rho = zeros(size(x,1),size(x,2),size(x,3),'single'); % < mean

    [s,rho,A,ll] = multicoil_infer(mag, ...
        'SensMaps',      s, ...
        'MeanImage',     rho, ...
        'Precision',     A, ...
        'RegCompFactor', prm, ...
        'VoxelSize',     vs, ...
        'IterMax',       5, ...
        'Tolerance',     1E-2, ...
        'SensOptim',     [true false], ...
        'CovOptim',      false, ...
        'Verbose',       2);

    s   = cat(5, s,   zeros(size(s),   'single'));
    rho = cat(5, rho, zeros(size(rho), 'single'));
else
    s   = zeros(size(x,1),size(x,2),size(x,3),size(x,4),2,'single');                     % < sensitivity
    rho = zeros(size(x,1),size(x,2),size(x,3),1,2,'single');
end

% -------------------------------------------------------------------------
% Main


[s,rho,A,ll] = multicoil_infer(x, ...
    'SensMaps',      s, ...
    'MeanImage',     rho, ...
    'Precision',     A, ...
    'RegCompFactor', prm, ...
    'VoxelSize',     vs, ...
    'Verbose',       2);