function [C,A] = multicoil_init_cov(x, gpu)
% Maximum likelihood covariance given a set of observed coil images, 
% log-sensitivity profiles and a mean image.
%
% FORMAT [C,A] = multicoil_init_cov(x, (gpu))
%
% x   - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
% gpu                                 - If 'gpu', compute on GPU
%
% C   -       Array [Nc Nc]           - Noise covariance matrix
% A   -       Array [Nc Nc]           - Noise precision matrix
%
% >> C = sum((b*rho-x)*(b*rho-x)')/(2J)     [where ' = conjugate transpose]
%
% Nc = number of coils
% Images can either be complex or have two real components that are then 
% assumed to be the real and imaginary parts.
% An output FileArray can be provided by using `rho` as an input. If not 
% provided, the output volume will have the same format as the input coil 
% volume.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 4
    gpu = '';
end
gpu_on = strcmpi(gpu, 'gpu');

% -------------------------------------------------------------------------
% Allocate output array
C = zeros(size(x,4), 'single');
if gpu_on
    C = gpuArray(C);
end
if gpu_on, loadarray = @loadarray_gpu;
else,      loadarray = @loadarray_cpu; end

N = size(x,4);
ndigits = ceil(log10(N))+1;

fprintf(['InitCov: %' num2str(ndigits) 's'], '');
for n=1:N
    
    fprintf([repmat('\b', [1 ndigits]) '%' num2str(ndigits) 'd'], n);
    
    % ---------------------------------------------------------------------
    % Load one coil image
    xn = loadarray(abs(x(:,:,:,1)), @single);
    xn = reshape(xn, [], 1);
    
    % ---------------------------------------------------------------------
    % Compute histogram
    xmax = max(xn,[],'omitnan');
    xmin = min(xn,[],'omitnan');
    cn = linspace(xmin,xmax,129);
    cn = (cn(2:end)+cn(1:end-1))/2;
    bwn = (xmax - xmin)/128;
    [xn,wn] = spm_imbasics('hist', double(xn), cn(:), 'KeepZero', false);
    clear c1
    
    % ---------------------------------------------------------------------
    % Fit Rician mixture
    
    [PI,NU,SIG] = spm_rice_mixture(double(wn), double(xn), 2);
    if NU(1) == 0 && NU(2) == 0
        [~,MU,A,PI] = spm_gmm(double(xn), 2, double(wn),...
            'GaussPrior', {[],10,[],10},'BinWidth',bwn);
        if MU(1) <= MU(2)
            C(n,n) = 1./A(1);
        else
            C(n,n) = 1./A(2);
        end
    else
        if NU(1) <= NU(2)
            C(n,n) = SIG(1)^2;
        else
            C(n,n) = SIG(2)^2;
        end
    end
    
end
fprintf('\n');

A = spm_matcomp('inv', C);