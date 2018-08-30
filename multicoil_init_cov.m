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

% -------------------------------------------------------------------------
% Allocate output array
C = zeros(size(x,4), 'single');
if strcmpi(gpu, 'gpu')
    C = gpuArray(C);
end

fprintf('InitCov:');
for n=1:size(x,4)
    
    fprintf(' %d', n);
    
    % ---------------------------------------------------------------------
    % Load one coil image
    if size(x, 5) == 2
        % Two real components
        x1 = reshape(single(x(:,:,:,n,:)), [], 2);
        x1 = sqrt(sum(x1.^2,2));
    else
        % One complex volume
        x1 = reshape(single(x(:,:,:,n)), [], 1);
        x1 = abs(x1);
    end
    
    % ---------------------------------------------------------------------
    % Compute histogram
    xmax = max(x1,[],'omitnan');
    xmin = min(x1,[],'omitnan');
    c1 = linspace(xmin,xmax,129);
    c1 = (c1(2:end)+c1(1:end-1))/2;
    bw1 = (xmax - xmin)/128;
    [x1,w1] = spm_imbasics('hist', double(x1), c1(:), 'KeepZero', false);
    clear c1
    
    % ---------------------------------------------------------------------
    % Fit Rician mixture
    
    [PI,NU,SIG] = spm_rice_mixture(double(w1), double(x1), 2);
    if NU(1) == 0 && NU(2) == 0
        [~,MU,A,PI] = spm_gmm(double(x1), 2, double(w1),...
            'GaussPrior', {[],10,[],10},'BinWidth',bw1);
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