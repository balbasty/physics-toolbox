function rho = multicoil_mean(varargin)
% Update mean given a set of observed coil images, 
% log-sensitivity profiles and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT rho = multicoil_mean(x, s, A, (rho), (prm), (vs))
%
% x   - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% A   -       Array [Nc Nc]           - Noise precision matrix
% rho - (File)Array [Nx Ny Nz  1 (2)] - Complex mean image [allocate]
% prm - 
% vs  - 
%
% If no input `rho` provided or `prm` == 0
%   => Maximum-Lkelihood estimate (closed-form)
% Else
%   => Maximum a poteriori estimate (Gauss-Newton)
%
% Nc = number of coils
% Images can either be complex or have two real components that are then 
% assumed to be the real and imaginary parts.
% An output FileArray can be provided by using `rho` as an input. If not 
% provided, the output volume will have the same format as the input coil 
% volume.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 5 || sum(varargin{5}) == 0
    if numel(varargin) > 4, varargin = varargin(1:4); end
    rho = multicoil_mean_ml(varargin{:});
else
    rho = multicoil_mean_map(varargin{:});
end