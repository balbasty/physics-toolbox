function dat = input(fname)
% Setup input data structure for Bloch-Siegert data.
% It is assuming a classic (single echo) acquisition.
%
% FORMAT dat = b1p.bs.io.input(fname)
% fname - Filename of h5 file
% dat   - Structure with fields:
%         . coils             - Reconstructed Bloch-Siegert data
%         . vs                - Voxel 
%         . pulse.voltage     - Voltage of the off-resonance pulse [in ?]
%         . pulse.voltage_ref - Reference voltage for 180deg/1s [in ?]
%         . pulse.frequency   - Offset frequency [in Hz]
%         . pulse.duration    - Off-resonance pulse duration [in s]
%         . pulse.shape       - Pulse shape (gauss|fermi)
%         . pulse.sign        - Sign of the pulse for each 'set' index
%         . pulse.factor      - Phase shift to B1 sqaured factor
%
% dat.coils is complex with shape:
%   [Nphase1 (Nphase2|Nslice) Nreadout Nchannel Nset]
if ~exist(fname, 'file')
    error('Input file does not seem to exist');
end
[~,~,ext] = fileparts(fname); 

switch lower(ext)
    case '.h5'
        dat = input_from_h5(fname);
    otherwise
        error('Input format not supported. File must be in ISMRMRD (h5) format.');
end
   

function dat = input_from_h5(fname)

dat = struct;

% -------------------------------------------------------------------------
% METADATA
% -------------------------------------------------------------------------
hdr = ismrmrd.xml(fname);
hdr = hdr.ismrmrdHeader;

% --- Bloch-Siegert parameters
% This section relies on user parameters specific to Nadege's sequence
fields                = [hdr.userParameters.userParameterDouble.name];
values                = [hdr.userParameters.userParameterDouble.value];
dat.pulse.voltage_ref = double(values(fields == "RefVoltage"));
dat.pulse.voltage     = double(values(fields == "BSPulseVoltage"));
dat.pulse.frequency   = double(values(fields == "BSPulse_OffsetFrequencyHz")); % Hz
fields                = [hdr.userParameters.userParameterLong.name];
values                = [hdr.userParameters.userParameterLong.value];
dat.pulse.duration    = double(values(fields == "BSPulseDuration")) * 1e-6; % sec
dat.pulse.sign        = [1 -1];
dat.pulse.shape       = 'fermi';
dat.pulse.factor      = b1p.bs.factor(dat.pulse.shape,     ...
                                      dat.pulse.frequency, ...
                                      dat.pulse.duration,  ...
                                      dat.pulse.voltage,   ...
                                      dat.pulse.voltage_ref);

% -------------------------------------------------------------------------
% DATA
% -------------------------------------------------------------------------

% --- Acquisition/Reconstruction metadata
acq_mtx = [hdr.encoding.encodedSpace.matrixSize.x ...
           hdr.encoding.encodedSpace.matrixSize.y ...
           hdr.encoding.encodedSpace.matrixSize.z];
acq_mtx = double(acq_mtx([2 3 1]));  % k1/(k2|sl)/rd
rec_mtx = [hdr.encoding.reconSpace.matrixSize.x ...
           hdr.encoding.reconSpace.matrixSize.y ...
           hdr.encoding.reconSpace.matrixSize.z];
rec_mtx = double(rec_mtx([2 3 1]));  % k1/(k2|sl)/rd
rec_fov = [hdr.encoding.reconSpace.fieldOfView_mm.x ...
           hdr.encoding.reconSpace.fieldOfView_mm.y ...
           hdr.encoding.reconSpace.fieldOfView_mm.z];
rec_fov = double(rec_fov([2 3 1]));  % k1/(k2|sl)/rd
dat.vs  = rec_fov./rec_mtx;

% --- Load k-space data (we don't know if 2D or 3D so we load k2 AND sl
dat.coils = ismrmrd.read(fname, 'contrast', 1, 'order', {'k1' 'k2' 'sl' 'rd' 'ch' 'st'});
dat.coils = utils.ifft(dat.coils, [1 2 4]); % IFFT along freq dimensions
dim       = size(dat.coils);
if dim(2) == 1
    dim(2) = [];
elseif dim(3) == 1
    dim(3) = [];
else
    error(['Strange... there seems to be both a second phase-encoding ' ...
           'direction AND a slice-encoding direction.']);
end
dat.coils = reshape(dat.coils, dim); % Remove unused dimension (k2 or sl)

% --- Crop if Recon < Acquisition
rec_off   = (acq_mtx-rec_mtx)/2;
dat.coils = dat.coils((rec_off(1)+1):(rec_off(1)+rec_mtx(1)), ...
                      (rec_off(2)+1):(rec_off(2)+rec_mtx(2)), ...
                      (rec_off(3)+1):(rec_off(3)+rec_mtx(3)), :,:,:,:);
                  
% TODO: Currently, crop assumes acquired frequencies are centred (no 
% partial Fourier) and even. It would be good to handle more general cases.
% Then, we should probably return a k-space undersampling mask as well.

% -------------------------------------------------------------------------
% NOISE ESTIMATE 
% -------------------------------------------------------------------------
dat.prec = zeros(size(dat.coils,4), size(dat.coils,5));
for c=1:size(dat.coils,4)
    for p=1:size(dat.coils, 5)
        dat.prec(c,p) = utils.noise_estimate(abs(dat.coils(:,:,:,c,p)));
        dat.prec(c,p) = 1./(dat.prec(c,p).^2);
    end
end