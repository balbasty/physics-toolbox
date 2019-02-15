function recon = multicoil_sense(varargin)

% -------------------------------------------------------------------------
% Helper functions
function ok = isarray(X)
    ok = islogical(X) || isnumeric(X) || isa(X, 'file_array');
end
function ok = isboolean(X)
    ok = (isnumeric(X) || islogical(X)) && isscalar(X);
end

% -------------------------------------------------------------------------
% Defaults
sens        = [];
A           = 1;
coilorder   = {'ch' 'rd' 'k1' 'k2' 'av' 'sl' 'ct' 'st' 'sg'};
sensorder   = {'k1' 'k2' 'rd' 'ch'};
af          = NaN;
coilcompact = true;
lat_acq     = NaN;
fov_acq     = NaN;
lat_recon   = NaN;
fov_recon   = NaN;
contrasts   = 1;
senscontrast= 1;
verbose     = 0;
parallel    = 0;

% -------------------------------------------------------------------------
% Parse input
p = inputParser;
p.FunctionName = 'multicoil_sense';
p.addRequired('CoilKSpace',                       @(X) ischar(X) || isarray(X));
p.addParameter('SensMaps',          sens,         @isarray);
p.addParameter('Precision',         A,            @isnumeric);
p.addParameter('CoilOrder',         coilorder,    @iscell);
p.addParameter('SensOrder',         sensorder,    @iscell);
p.addParameter('CoilCompact',       coilcompact,  @isboolean);
p.addParameter('Acceleration',      af,           @isnumeric);
p.addParameter('AcquisitionMatrix', lat_acq,      @isnumeric);
p.addParameter('AcquisitionFOV',    fov_acq,      @isnumeric);
p.addParameter('ReconMatrix',       lat_recon,    @isnumeric);
p.addParameter('ReconFOV',          fov_recon,    @isnumeric);
p.addParameter('Contrast',          contrasts,    @isnumeric);
p.addParameter('SensContrast',      senscontrast, @isnumeric);
p.addParameter('Parallel',          parallel,     @isboolean);
p.addParameter('Verbose',           verbose,      @isboolean);
p.parse(varargin{:});
rdata        = p.Results.CoilKSpace;
sens         = p.Results.SensMaps;
af           = p.Results.Acceleration;
A            = p.Results.Precision;
coilorder    = p.Results.CoilOrder;
sensorder    = p.Results.SensOrder;
coilcompact  = p.Results.CoilCompact;
lat_acq      = p.Results.AcquisitionMatrix;
fov_acq      = p.Results.AcquisitionFOV;
lat_recon    = p.Results.ReconMatrix;
fov_recon    = p.Results.ReconFOV;
contrasts    = p.Results.Contrast;
senscontrast = p.Results.SensContrast(:)';
parallel     = p.Results.Parallel;
verbose      = p.Results.Verbose;


% Parallel processing
if isnan(parallel)
    parallel = 0;
end
if ~isfinite(parallel)
    parallel = parcluster('local');
    parallel = parallel.NumWorkers;
end

% -------------------------------------------------------------------------
% Estimate senstivity
if isempty(sens)
    if ischar(rdata)
        acdata = ismrmrd_read(rdata, ...
                              'contrast', senscontrast, ...
                              'subpart',  'autocalib');
        % order: [ch rd k1 k2]
        acdata = fftshift(ifft(ifftshift(fftshift(ifft(ifftshift(fftshift(ifft(ifftshift(acdata,2),[],2),2),3),[],3),3),4),[],4),4);
        acdata = permute(acdata, [3 4 2 1]); % < order: [k1 k2 rd ch]
        acdim  = size(acdata);
        sens   = zeros([acdim 2],        'like', acdata);
        acmean = zeros([acdim(1:3) 1 2], 'like', acdata);
        [C,A] = multicoil_init_cov(acdata);
        [sens,acmean,A,ll] = multicoil_infer(acdata, ...
            'SensMaps',      sens, ...
            'MeanImage',     acmean, ...
            'Precision',     A, ...
            'RegCompFactor', 1E7, ...
            'VoxelSize',     [4.48 5.6 6.4], ...
            'Verbose',       1);
    else
        error('Calibration data needed to estimate sensitivity maps')
    end
end

% -------------------------------------------------------------------------
% Convert log-sensitivity to frequency (Neumann conditions -> dct)
sens    = permute(sens, [4 1 2 3]); % < order: [ch k1 k2 rd]
sens    = dct(dct(dct(sens, [], 2), [], 3), [], 4);
dimsens = size(sens);
dimsens = dimsens(2:4);

for contrast=contrasts

% -------------------------------------------------------------------------
% Read and reorganize acquired data
if ischar(rdata)
    fname = rdata;
    rdata = ismrmrd_read(fname, ...
                         'contrast', contrast, ...
                         'subpart',  'cartesian', ...
                         'compact',  true);
    % order: [ch rd k1 k2]
end
% recon_lat = reconstruction matrix
recon_lat    = size(rdata);
recon_lat(1) = 1;
recon_lat(3) = recon_lat(3) * af(1);
recon_lat(4) = recon_lat(4) * af(2);
recon_lat    = recon_lat([3 4 2]);
rdata        = permute(rdata, [1 3 4 2]); % < order: [ch k1 k2 rd]

% -------------------------------------------------------------------------
% Prepare a bit of stuff for sense

% > Sampling scheme (k1/k2 are transposed compared to rdata order)
msk = zeros(recon_lat(2),recon_lat(1),'logical'); % [k2 k1]
msk(1:af(2):end,1:af(1):end) = 1;

% > Aliasing locations and coefficients
msksp   = ifft2(ifftshift(msk));
[k2,k1] = find(abs(msksp) > eps('single'));
weights = double(msksp(abs(msksp) > eps('single')));

% > ifft in readout direction (done once and for all)
rdata = fftshift(ifft(ifftshift(rdata,4), [], 4), 4);
% > idct in readout direction (done once and for all)

sens = idct(sens, recon_lat(3), 4);
sens = sens * sqrt(recon_lat(3)) / sqrt(dimsens(3));

% % > debug: test on a few slices
% rdata = rdata(:,:,:,76:95);
% sens  = sens(:,:,:,76:95,:);
% sensmsk  = sensmsk(:,:,:,76:95);
% ssq  = ssq(:,:,:,76:95);
% recon_lat(3) = 20;

% -------------------------------------------------------------------------
% Reconstruction
recon = zeros(recon_lat(1)*recon_lat(2),recon_lat(3), 'like', single(1i));
ncoil = size(rdata,1);
parfor(z=1:recon_lat(3), parallel)
    % Acquired: read one slice (and expand to full matrix)
    xz = zeros(ncoil,recon_lat(1),recon_lat(2),'like', double(1i));
    xz(:,1:af(1):end,1:af(2):end) = rdata(:,:,:,z);
    % Inverse Fourier to generate aliased image
    xz = fftshift(ifft(ifftshift(fftshift(ifft(ifftshift(xz,2),[],2),2),3),[],3),3);
    % Sensitivity: read one slice 
    sz = double(sens(:,:,:,z));
    % Frequency -> Image
    sz = idct(idct(sz,recon_lat(1),2),recon_lat(2),3);
    sz = sz * sqrt(recon_lat(1)) / sqrt(dimsens(1)) ...
            * sqrt(recon_lat(2)) / sqrt(dimsens(2));
    % Exponentiate
    sz = exp(sz);
    % Whiten
    xz = reshape(xz, ncoil, []);
    sz = reshape(sz, ncoil, []);
    xz = sqrt(A) * xz;
    sz = sqrt(A) * sz;
    xz = reshape(xz, [ncoil recon_lat(1:2)]);
    sz = reshape(sz, [ncoil recon_lat(1:2)]);
    % SENSE inversion
    recon(:,z) = SENSEUnfold_v5(xz, sz, [k1 k2], [af(2) af(1)], false, weights);
end
recon = reshape(recon,recon_lat);
recon = permute(recon, [3 1 2]); % < order: [rd k1 k2]

end % < contrast loop

end % < function