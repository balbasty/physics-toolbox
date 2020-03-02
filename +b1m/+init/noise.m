function [A,N] = noise(coils, mask)
% Estimate the noise precision of each coil by fitting a Rice or Gaussian
% mixture model to their magnitude image.
%
% FORMAT [prec,n] = b1m.init.noise(coils, (mask))
%
% coils - (File)Array [Nx Ny Nz Nc] - Complex coil images
% mask  - (File)Array [Nx Ny]       - Sampling mask
%
% prec  -       Array [Nc Nc]       - Diagonal noise precision matrix
% n     -                           - Number of points in image domain
%
% This function returns noise precision in the image domain of the
% fully-sampled centre of k-space. The precision of the noise in k-space is
% obtained by: prec/n
%
% Nc = number of coils
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 2, mask = logical([]); end

% -------------------------------------------------------------------------
% Get autocalibration region
if numel(mask) > 1
    coils = utils.ifft(utils.acsub(coils, mask), [1 2 3]);
end
N = size(coils,1)*size(coils,2)*size(coils,3);

% -------------------------------------------------------------------------
% Allocate output array
C = zeros(size(coils,4), 'single');

Nc = size(coils,4);
ndigits = ceil(log10(Nc))+1;

fprintf(['InitCov: %' num2str(ndigits) 's'], '');
for n=1:Nc
    
    fprintf([repmat('\b', [1 ndigits]) '%' num2str(ndigits) 'd'], n);
    
    % ---------------------------------------------------------------------
    % Load one coil image
    xn = single(abs(coils(:,:,:,n)));
    xn = reshape(xn, [], 1);
    
    % ---------------------------------------------------------------------
    % Compute histogram
    xmax = max(xn,[],'omitnan');
    xmin = min(xn,[],'omitnan');
    cn = linspace(xmin,xmax,129);
    cn = (cn(2:end)+cn(1:end-1))/2;
    bwn = (xmax - xmin)/128;
    [xn,wn] = utils.histN(double(xn), cn(:), 'KeepZero', false);
    clear c1
    
    % ---------------------------------------------------------------------
    % Fit Rician mixture
    [~,~,SIG] = spm_rice_mixture(double(wn), double(xn), 2);
    SIG = min(SIG);
    if SIG == 0
        % If failed, fit Gaussian mixture
        [~,~,A] = utils.gmm.fit(double(xn), 2, double(wn),...
            'GaussPrior', {[],10,[],10},'BinWidth',bwn);
        C(n,n) = 1./max(A);
    else
        C(n,n) = SIG.^2;
    end
    
end
fprintf('\n');

A = utils.matcomp('inv', C);