% -------------------------------------------------------------------------
% Parameters
vs    = [8 8 8];    % > Voxel size
itmax = 1000;       % > Maximum number of optim iterations
prm   = [0 0 1E7];  % > Sensitivity fields regularisation 
                    %   ['abs value' 'membrane E' 'bending E']

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
% N = 8;
x = x(:,:,:,idx(1:N),:);

% -------------------------------------------------------------------------
% Allocate unknown arrays
s   = zeros(size(x),'single');                                   % < sensitivity
rho = zeros(size(x,1),size(x,2),size(x,3),1,size(x,5),'single'); % < mean

% -------------------------------------------------------------------------
% Noise estimates
% > We start with some sort of initial estimate based on a very
%   simple independant Gaussian (mixture) fitting
[C,A] = multicoil_init_cov(x); % < C = covariance / A = inv(C) = precision

% -------------------------------------------------------------------------
% Main
[s,rho,A,ll] = multicoil_infer(x, s, rho, A, prm, vs, itmax, 2);