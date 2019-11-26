function [sens,meanim,prec] = multiscale(varargin)
% Compute mode estimates (ML, MAP) of the parameters of a probabilistic 
% model of complex multicoil MR images, using a multiscale (multigrid) 
% approach.
%
% FORMAT [sens,mean,prec,ll] = b1m.multiscale(coils, ...)
%
% REQUIRED
% --------
% coils - (File)Array [Nx Ny Nz Nch Nct] - Complex coil k-spaces
%
% KEYWORDS
% --------
% Precision         - Noise precision matrix                    [NaN=estimate from fully-sampled]
% SensRegFactor     - Regularisation factor                     [1E2]
% Relax             - Relaxation iterations                     [5]
% Recon             - Reconstruction matrix                     [NaN=same as input]
% VoxelSize         - Voxel size                                [1]
% Verbose           - -1  = quiet
%                     [0] = moderate
%                      1  = verbose
%                      2  = plot mean image and log-likelihood
%                      3  = plot coil-wise fit (slow + disables parallel processing)
%
% ADVANCED
% --------
% Parallel          - Number of parallel workers                [Inf=all]
%
% OUTPUT
% ------
% sens - (Log)-Sensitivity maps - (File)Array [Nx Ny Nz Nc]
% mean - Mean image             - (File)Array [Nx Ny Nz Nc]
% prec - Noise precision        -       Array [Nc Nc]
%
% An output FileArray can be provided by using `MeanImage` as an input. If  
% not provided, the output volume will have the same format as the input  
% coil volume.
%
% Nch = number of coils
% Nx  = Phase encode 1 
% Ny  = Phase encode 2 /or/ Slice
% Nz  = Frequency readout
% Nct = number of contrasts
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

start = tic;

% =========================================================================
%
%                       PARSE AND PROCESS ARGUMENTS
%
% =========================================================================

% -------------------------------------------------------------------------
% Parse input
% -------------------------------------------------------------------------
p  = inputParser;
p.FunctionName = 'b1m.multiscale';
p.addRequired('CoilImages',                      @utils.isarray);
p.addParameter('Precision',         NaN,         @isnumeric);
p.addParameter('SensRegFactor',     1E2,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('VoxelSize',         [1 1 1],     @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('Relax',             5,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('Recon',             NaN,         @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('Parallel',          Inf,         @utils.isboolean);
p.addParameter('Verbose',           0,           @utils.isboolean);
p.parse(varargin{:});
coils           = p.Results.CoilImages;
prec            = p.Results.Precision;
relax           = p.Results.Relax;
sensfactor      = p.Results.SensRegFactor;
vs              = p.Results.VoxelSize;
verbose         = p.Results.Verbose;
Nw              = p.Results.Parallel;
recon           = p.Results.Recon;

% -------------------------------------------------------------------------
% Problem size
% -------------------------------------------------------------------------
Nx  = size(coils,1);
Ny  = size(coils,2);
Nz  = size(coils,3);
Nc  = size(coils,4);
Nct = size(coils,5);
lat = [Nx Ny Nz];

scales = {[Nx Ny Nz]};
while all(scales{end} > 16)
    scales = [scales {ceil(scales{end}/4)*2}];
end

recon = utils.pad(recon(:)', [0 3-numel(recon)], 'replicate', 'post');
recon = recon(1:3);
recon(~isfinite(recon)) = lat(~isfinite(recon));

% -------------------------------------------------------------------------
% Initialise scales
% -------------------------------------------------------------------------
coils  = [{coils} cell(1,numel(scales)-1)];
sens   = cell(1,numel(scales));
[sens{:}] = deal([]);
meanim = cell(1,numel(scales));
[meanim{:}] = deal([]);
mask = cell(1,numel(scales));
[mask{:}] = deal([]);
for i=1:numel(scales)
    if i > 1
        coils{i} = utils.ksub(coils{i-1}, scales{i});
    end
    mask{i}  = sum(sum(sum(abs(coils{i}), 3), 4), 5) > single(eps);
end

% -------------------------------------------------------------------------
% Compute precision
% -------------------------------------------------------------------------
if any(~isfinite(prec(:)))
    ac = utils.ifft(utils.acsub(coils{1}, mask{1}), [1 2 3]);
    prec = b1m.init.noise(ac);
    prec = prec/(size(ac,1)*size(ac,2)*size(ac,3));
end

% -------------------------------------------------------------------------
% Fourier transform
% -------------------------------------------------------------------------
for i=1:numel(scales)
    coils{i}  = utils.ifft(coils{i}, [1 2 3]);
    scales{i} = ceil(scales{i}.*recon./scales{1});
    coils{i}  = utils.ksub(coils{i}, scales{i});
end

% -------------------------------------------------------------------------
% Multigrid
% -------------------------------------------------------------------------
i        = numel(scales);
while i > 0
    
    fprintf('. Level %d: [%d %d %d]\n', ...
        i, scales{i}(1), scales{i}(2), scales{i}(3));
    
        
    if numel(scales) > i
        meanim{i} = complex(spm_diffeo('resize', real(single(meanim{i+1})), scales{i}), ...
                            spm_diffeo('resize', imag(single(meanim{i+1})), scales{i}));
%         sens{i}   = b1m.upsample(sens{i+1}, scales{i});
        sens{i}   = complex(spm_diffeo('resize', real(single(sens{i+1})), scales{i}), ...
                            spm_diffeo('resize', imag(single(sens{i+1})), scales{i}));
        meanim{i} = meanim{i} / prod(scales{i}) * prod(scales{i+1});
    end
    
    [sens{i},meanim{i},prec] = b1m.fit(coils{i}, ...
        'SamplingMask',     mask{i}, ...
        'SensMaps',         sens{i}, ...
        'MeanImage',        meanim{i}, ...
        'VoxelSize',        vs*2^(i-1), ...
        'Parallel',         Nw, ...
        'SensRegFactor',    sensfactor*(prod(scales{1})/prod(scales{i})).^2, ...
        'MeanRegfactor',    0, ... 
        'SensLog',          true, ...
        'IterMax',          relax, ...
        'IterMin',          relax, ...
        'IterDec',          0, ...
        'Precision',        prec, ...
        'Verbose',          verbose);
        
    
    % Move up
    fprintf('Move up!\n');
    scales = scales(1:i);
    sens   = sens(1:i);
    meanim = meanim(1:i);
    coils  = coils(1:i);
    mask   = mask(1:i);
    
    i = i-1;
    
end

meanim = meanim{1};
sens   = sens{1};

stop = toc(start);
fprintf('Total time: %s\n', utils.sec2ydhms(stop));