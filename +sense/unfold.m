function recon = unfold(varargin)
% Unfold a subsampled acquisition using pre-computed sensitivity fields.
%
% FORMAT recon = SENSE.unfold(coils, sens, ...)
%
% REQUIRED
% --------
% coils - (File)Array [Nc Nx Ny Nz Nct] - Undersampled k-space data
% sens  - (File)Array [Nc Nx Ny Nz Nct] - Sensitivity fields
%
% Nc  = Coils
% Nx  = Phase-encoding direction 1
% Ny  = Phase-encoding direction 2 /or/ slice
% Nz  = Frequency readout
% Nct = Contrasts
%
% KEYWORD
% -------
% CoilCompact  - Is k-space stored compactely?              [true]
% CoilSpace    - Frequency (=0) or image (=1) space of      [0]
%                each dimension
% SensLog      - Sensitivities in log-space?                [true]
% SensPInv     - Are sensitivities already pseudo-inverted? [false]
% ReconMatrix  - Reconstruction matrix                      [Nan=auto]
% Acceleration - Acceleration factor [k1 k2]                [Nan=auto]
% Precision    - Noise precision matrix                     [1=identity]
% Contrast     - Contrast to unfold                         [Inf=all]
% Parallel     - Number of parallel workers                 [Inf]
% Verbose      - Verbosity level                            [0]
% 
% OUTPUT
% ------
% recon - (File)Array [Nx Ny Nz Nct] - Reconstructed volume
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% =========================================================================
%
%                       PARSE AND PROCESS ARGUMENTS
%
% =========================================================================

% -------------------------------------------------------------------------
% Defaults
coil_compact = true;
coil_space   = 0;
sens_log     = true;
sens_pinv    = false;
recon_lat    = NaN;
af           = NaN;
precision    = 1;
contrasts    = 1;
verbose      = 0;
parallel     = Inf;

% -------------------------------------------------------------------------
% Parse input
p = inputParser;
p.FunctionName = 'sense.unfold';
p.addRequired('CoilKSpace',                  @utils.isarray);
p.addRequired('SensMaps',                    @utils.isarray);
p.addParameter('CoilCompact',  coil_compact, @utils.isboolean);
p.addParameter('CoilSpace',    coil_space,   @isnumeric);
p.addParameter('SensLog',      sens_log,     @utils.isboolean);
p.addParameter('SensPInv',     sens_pinv,    @utils.isboolean);
p.addParameter('ReconMatrix',  recon_lat,    @isnumeric);
p.addParameter('Acceleration', af,           @isnumeric);
p.addParameter('Precision',    precision,    @isnumeric);
p.addParameter('Contrast',     contrasts,    @isnumeric);
p.addParameter('Parallel',     parallel,     @utils.isboolean);
p.addParameter('Verbose',      verbose,      @utils.isboolean);
p.parse(varargin{:});
rdata        = p.Results.CoilKSpace;
sens         = p.Results.SensMaps;
coil_compact = p.Results.CoilCompact;
coil_space   = p.Results.CoilSpace;
sens_log     = p.Results.SensLog;
sens_pinv    = p.Results.SensPInv;
recon_lat    = p.Results.ReconMatrix;
af           = p.Results.Acceleration;
precision    = p.Results.Precision;
contrasts    = p.Results.Contrast;
verbose      = p.Results.Verbose;
parallel     = p.Results.Parallel;


% -------------------------------------------------------------------------
% Preprocess arguments

% Sanity check dimensions
if size(rdata,1) ~= size(sens,1)
    error(['Number of coils is different between sensitivities and ' ...
           'k-space data (%d ~= %d)'], size(sens,1), size(rdata,1));
end
Nc = size(rdata, 1); % Number of coils

% Pad coil space
if isempty(coil_space), coil_space = 0; end
coil_space = padarray(coil_space(:)', [0 max(0, 3-numel(coil_space))], ...
                      'replicate', 'post');

% Precision
if isempty(precision), precision = 1; end
if size(precision,1) == 1 || size(precision,2) == 1
    precision = padarray(precision(:)', [0 max(Nc,Nc-numel(precision))], ...
                         'replicate', 'post');
    precision = diag(precision);
end

% Parallel processing
if isnan(parallel)
    parallel = 0;
end
if ~isfinite(parallel)
    parallel = parcluster('local');
    parallel = parallel.NumWorkers;
end

% Sanity check on SensLog/SensPInv
if sens_log && sens_pinv
    warning(['Options Senslog and SensPInv are both set to true, '    ...
             'which is contradictory. I will assume SensPInv = true ' ...
             'and SensLog = false']);
    sens_log = false;
end

% Reconstruction matrix
recon_lat = padarray(recon_lat(:)', [0 max(0, 3-numel(recon_lat))], ...
                     'replicate', 'post');
if ~isfinite(recon_lat(1)) % k1
    if ~coil_compact
        recon_lat(1) = size(rdata,2);
    else
        if isfinite(af(1))
            recon_lat(1) = size(rdata,2)*af(1);
            if mod(recon_lat(1),2) == 1
                recon_lat(1) = recon_lat(1) + 1;
            end
        else
            error(['Cannot determine the reconstruction matrix based ' ...
                   'on the provided inputs.']);
        end
    end
end
if ~isfinite(recon_lat(2)) % k2
    if ~coil_compact
        recon_lat(2) = size(rdata,3);
    else
        if isfinite(af(2))
            recon_lat(2) = size(rdata,3)*af(2);
            if mod(recon_lat(2),2) == 1
                recon_lat(2) = recon_lat(2) + 1;
            end
        else
            error(['Cannot determine the reconstruction matrix based ' ...
                   'on the provided inputs.']);
        end
    end
end
if ~isfinite(recon_lat(3)) % rd
    recon_lat(3) = size(rdata,4);
end

% -------------------------------------------------------------------------
% Sampling scheme (k1/k2 are transposed compared to internal order)
msk = zeros(recon_lat(2),recon_lat(1),'logical'); % [k2 k1]
msk(1:af(2):end,1:af(1):end) = 1;

% -------------------------------------------------------------------------
% Aliasing locations and coefficients
msksp   = ifft2(ifftshift(msk));
[k2,k1] = find(abs(msksp) > eps('single'));
weights = double(msksp(abs(msksp) > eps('single')));

% -------------------------------------------------------------------------
% idct in readout direction (done once and for all)
sens = b1m.upsample(single(sens()), [NaN NaN NaN recon_lat(3)]);

% =========================================================================
%
%                          LOOP OVER CONTRASTS
%
% =========================================================================

for contrast=contrasts

if verbose > 0
    fprintf('Unfolding contrast %d...\n', contrast);
end

% -------------------------------------------------------------------------
% Load one contrast
rdata1 = single(rdata(:,:,:,:,contrast));

% -------------------------------------------------------------------------
% ifft in readout direction (done once and for all)
if ~coil_space(3)
    rdata1 = utils.ifft(rdata1,4);
end

% -------------------------------------------------------------------------
% Reconstruction
recon = zeros(recon_lat(1)*recon_lat(2),recon_lat(3), 'like', single(1i));
parfor(z=1:recon_lat(3), parallel)
    % K-space
    if coil_compact
        xz = zeros([Nc recon_lat(1:2)],'like', double(1i));      % Allocate
        xz(:,1:af(1):end,1:af(2):end) = rdata1(:,:,:,z);
    else
        xz = rdata1(:,:,:,z);
    end
    if ~coil_space(1), xz = utils.ifft(xz, 2); end
    if ~coil_space(2), xz = utils.ifft(xz, 3); end
    % Sensitivity: 
    sz = double(sens(:,:,:,z));                            % Read one slice
    sz = b1m.upsample(sz, [NaN recon_lat(1:2)]);       % Frequency -> Image
    if sens_log, sz = exp(sz); end                           % Exponentiate
    % Whiten
    xz = reshape(xz, Nc, []);
    sz = reshape(sz, Nc, []);
    xz = sqrt(precision) * xz;
    sz = sqrt(precision) * sz;
    xz = reshape(xz, [Nc recon_lat(1:2)]);
    sz = reshape(sz, [Nc recon_lat(1:2)]);
    % SENSE inversion
    recon(:,z) = SENSEUnfold_v5(xz, sz, [k1 k2], af(2:-1:1), sens_pinv, weights);
end
recon = reshape(recon,recon_lat);


end % < contrast loop

end % < function