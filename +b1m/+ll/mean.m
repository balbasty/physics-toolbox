function llp = mean(varargin)
% Compute prior log-likelihood of the mean
%
% FORMAT ll = b1m.ll.mean(mean, ...)
%
% REQUIRED
% --------
% mean - (File)Array [Nx Ny Nz 1 Nct] - Complex mean image
%
% KEYWORD
% -------
% RegFactor - Scalar       - Regularisation factor      [1]
% VoxelSize - Array [1 3]  - Voxel size                 [1]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Parse input
p = inputParser;
p.FunctionName = 'b1m.ll.mean';
p.addRequired('MeanImage',                   @utils.isarray);
p.addParameter('RegFactor',     1,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('VoxelSize',     1,           @(X) isnumeric(X) && isrow(X) && numel(X) <= 3);
p.parse(varargin{:});
meanim  = p.Results.MeanImage;
reg     = p.Results.RegFactor;
vs      = p.Results.VoxelSize;

% -------------------------------------------------------------------------
% Regularisation structure
if reg > 0
    prm = [0 1 0];              % Membrane energy
    spm_field('boundary', 1);   % Neumann boundary conditions
else
    llp = 0;
    return
end

% -------------------------------------------------------------------------
% Compute log-likelihood (prior)
meanim = single(gather(meanim));
meanim = cat(4, real(meanim), imag(meanim));
llp    = spm_field('vel2mom', meanim, [vs reg*prm]);
llp    = -0.5 * double(meanim(:))' * double(llp(:));