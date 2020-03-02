function ll = sensitivity(varargin)
% Compute prior log-likelihood of the sensitivities
%
% FORMAT lls = b1m.ll.sensitivity(sens, ...)
%
% REQUIRED
% --------
% sens - (File)Array [Nx Ny Nz Nc] - Complex log-sensitivity profiles          
%
% KEYWORD
% -------
% RegFactor - Array [Nc 2] - Regularisation factor      [1]
% VoxelSize - Array [1 3]  - Voxel size                 [1]
%
% OUTPUT
% ------
% lls - Array [Nc 2] - Non-constant term of the log-likelihood, 
%                      per coil and part
%
% NOTE: RegFactor/VoxelSize are implicitely expanded along singleton dim.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

Nc = size(varargin{1}, 4);

% -------------------------------------------------------------------------
% Parse input
p = inputParser;
p.FunctionName = 'b1m.ll.sensitivity';
p.addRequired('SensMaps',                    @utils.isarray);
p.addParameter('RegFactor',     1,           @(X) isnumeric(X) && size(X,1) <= Nc && size(X,2) <= 2);
p.addParameter('VoxelSize',     1,           @(X) isnumeric(X) && isrow(X) && numel(X) <= 3);
p.parse(varargin{:});
sens    = p.Results.SensMaps;
reg     = p.Results.RegFactor;
vs      = p.Results.VoxelSize;

% -------------------------------------------------------------------------
% Pad parameters
reg = utils.pad(double(reg), [Nc-size(reg,1) 2-size(reg,2)], 'replicate', 'post');
vs  = utils.pad(double(vs), [0 3-numel(vs)], 'replicate', 'post');

% -------------------------------------------------------------------------
% Regularisation structure
regstruct = [0 0 1];       % Bending energy
spm_field('boundary', 1);  % Neumann boundary conditions
prm = [vs regstruct];

% -------------------------------------------------------------------------
% Compute log-likelihood (prior)
ll = zeros(Nc,2);
for n=1:Nc
    reg1    = reg(n,:);
    s1      = single(gather(sens(:,:,:,n)));
    sr      = real(s1);
    si      = imag(s1);
    s1      = [];
    llpr    = spm_field('vel2mom', sr, prm, reg1(1));
    llpr    = -0.5 * double(llpr(:))' * double(sr(:));
    sr      = [];
    llpi    = spm_field('vel2mom', si, prm, reg1(2));
    llpi    = -0.5 * double(llpi(:))' * double(si(:));
    si      = [];
    ll(n,:) = [llpr llpi];
end