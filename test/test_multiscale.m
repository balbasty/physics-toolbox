% ------------
% Load dataset
dataset = load_data;
c = 1;

% -------------
% Read filename
input_folder  = dataset(1).folder;
subfolder     = dataset(1).subfolder;
fname_array   = fullfile(input_folder, subfolder, dataset(1).T1w(1).array);

% -----------
% Read header
header        = ismrmrd.xml(fname_array);
header        = header.ismrmrdHeader;
recon_fov_rd  = header.encoding.reconSpace.fieldOfView_mm.x;
recon_fov_k1  = header.encoding.reconSpace.fieldOfView_mm.y;
recon_fov_k2  = header.encoding.reconSpace.fieldOfView_mm.z;
recon_fov     = [recon_fov_rd recon_fov_k1 recon_fov_k2];
recon_mtx_rd  = header.encoding.reconSpace.matrixSize.x;
recon_mtx_k1  = header.encoding.reconSpace.matrixSize.y;
recon_mtx_k2  = header.encoding.reconSpace.matrixSize.z;
recon_mtx     = [recon_mtx_rd recon_mtx_k1 recon_mtx_k2];
recon_vs      = recon_fov ./ recon_mtx;

% ---------
% Read data
kdata         = single(ismrmrd.read(fname_array, 'contrast', 1));
kdata         = kdata * 1000 / max(abs(kdata(:)));

% ---------------
% Reorganise data
kdata         = permute(kdata, [3 4 2 1 5]);
vs            = recon_vs([2 3 1]);

% -----
% Recon
[sens,meanim,prec] = b1m.multiscale(kdata, ...
    'SensRegFactor', 1E5, 'VoxelSize', vs, 'Verbose', 2);

% -----------------------------------------------
% Remove non-image AC lines + undersample readout
kdata         = utils.ifft(kdata, 3);
kdata         = utils.ksub(kdata, [NaN NaN recon_mtx_rd]);
kdata(2:2:end,:,:,:) = 0;
kdata(:,2:2:end,:,:) = 0;

% ---------------
% SENSE unfolding
recon = sense.unfold(permute(kdata, [4 1 2 3]), permute(sens, [4 1 2 3]), ...
    'CoilCompact',  false,    ...
    'CoilSpace',    [0 0 1], ...
    'Acceleration', [2 2], ...
    'Precision',    prec, ...
    'Parallel',     Inf);

ssq = recon .* sqrt(sum(sens.^2,4));