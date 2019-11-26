function [coils, prec, mask] = prepare(varargin)
% Prepare calibration data to be used with b1m.fit
%
% FORMAT [coils, prec, mask] = b1m.prepare(kdata, ...)
%
% REQUIRED
% --------
% kdata - [ch rd k1 k2 ct] - Calibration k-space data
%
% KEYWORD
% -------
% AcqFOV      - [rd k1 k2] - Acquisition field of view (mm)   [NaN=1mm iso]
% ReconFOV    - [rd k1 k2] - Reconstruction field of view     [NaN=AcqFOV]
% ReconMatrix - [rd k1 k2] - Calibration lattice              [NaN=same]
% FFT                      - Apply FFT?                       [true=all]
%
% OUTPUT
% ------
% coils - [k1 k2 rd ch ct] - Calibration image data
% prec  - [Nc Nc]          - Noise precision matrix
% mask  - [k1 k2]          - Sampling pattern
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Parse input
% -------------------------------------------------------------------------
p  = inputParser;
p.FunctionName = 'b1m.prepare';
p.addRequired('kdata',                  @utils.isarray);
p.addParameter('AcqFOV',      NaN,      @isnumeric);
p.addParameter('ReconFOV',    NaN,      @isnumeric);
p.addParameter('ReconMatrix', NaN,      @isnumeric);
p.addParameter('FFT',         true,     @(X) isnumeric(X) || islogical(X));
p.parse(varargin{:});
coils     = p.Results.kdata;
acq_fov   = p.Results.AcqFOV;
recon_fov = p.Results.ReconFOV;
recon_lat = p.Results.ReconMatrix;
do_fft    = p.Results.FFT;

% Pad arguments
acq_fov   = utils.pad(acq_fov(:)', [0 max(0, 3-numel(acq_fov))], ...
                     'replicate', 'post');
recon_lat = utils.pad(recon_lat(:)', [0 max(0, 3-numel(recon_lat))], ...
                     'replicate', 'post');
recon_fov = utils.pad(recon_fov(:)', [0 max(0, 3-numel(recon_fov))], ...
                     'replicate', 'post');
do_fft    = utils.pad(do_fft(:)', [0 max(0, 3-numel(do_fft))], ...
                     'replicate', 'post');

% Get acquisition lattice
Nc = size(coils,1); % Number of coils
Nr = size(coils,2); % Number of frequency encoding
N1 = size(coils,3); % Number of phase encoding 1
N2 = size(coils,4); % Number of phase encoding 2
Ne = size(coils,5); % Number of echos/contrasts
acq_lat = [Nr N1 N2];

% Default values
for i=1:3
    if ~isfinite(acq_fov(i))
        acq_fov(i) = acq_lat(i);
    end
end
for i=1:3
    if ~isfinite(recon_fov(i))
        recon_fov(i) = acq_fov(i);
    end
end
for i=1:3
    if ~isfinite(recon_lat(i))
        recon_lat(i) = ceil(acq_lat(i)*recon_fov(i)/acq_fov(i));
    end
end

% Estimate noise precision
if nargout >= 2
    [~,prec] = b1m.init.noise(permute(utils.ifft(coils, [2 3 4]), [3 4 2 1 5]));
    prec     = prec/prod(acq_lat);
end

% Prepare subsampling mask
if nargout >= 3
    mask = ones(acq_lat(2:3), 'logical');
end

fov_factor = acq_fov./recon_fov;

% % Subsample dimensions where recon_fov < acq_fov
% for i=1:3
%     if (recon_fov(i) < acq_fov(i)) && ~do_fft(i)
%         step   = ceil(acq_fov(i)/recon_fov(i));
%         S      = struct;
%         S.type = '()';
%         centre = floor(acq_lat(i)/2) + 1;
%         S.subs = repmat({':'}, [1 5]);
%         idx    = [fliplr(centre:-step:1) (centre+step):step:acq_lat(i)];
%         S.subs{i+1} = idx;
%         coils  = subsref(coils, S);
%         if nargout >= 3 && i ~= 1
%             S.subs = S.subs([3 4]);
%             mask = subsref(mask, S);
%         end
%     end
% end

% Update acquisition lattice
Nr = size(coils,2);
N1 = size(coils,3);
N2 = size(coils,4);
acq_lat = [Nr N1 N2];

recon_lat = recon_lat .* fov_factor;

% Crop dimensions where recon_lat < acq_lat
for i=1:3
    if recon_lat(i) < acq_lat(i)
        crop   = ceil(acq_lat(i) - recon_lat(i));
        S      = struct;
        S.type = '()';
        if mod(acq_lat(i),2)
            % odd lattice
            if mod(crop,2)
                % odd crop
                idx = ceil(crop/2):(acq_lat(i)-ceil(crop/2));
            else
                % even crop
                idx = (crop/2+1):(acq_lat(i)-(crop/2));
            end
        else
            % even lattice
            if mod(crop,2)
                % odd crop
                idx = ceil(crop/2+1):(acq_lat(i)-floor(crop/2));
            else
                % even crop
                idx = (crop/2+1):(acq_lat(i)-(crop/2));
            end
        end
        S.subs = repmat({':'}, [1 5]);
        S.subs{i+1} = idx;
        coils  = subsref(coils, S);
        if nargout >= 3 && i ~= 1
            S.subs = S.subs([3 4]);
            mask = subsref(mask, S);
        end
    end
end
% Update acquisition lattice
Nr = size(coils,2);
N1 = size(coils,3);
N2 = size(coils,4);
acq_lat = [Nr N1 N2];

% Pad dimensions where recon_lat > acq_lat
for i=1:3
    if recon_lat(i) > acq_lat(i)
        diff     = ceil(recon_lat(i) - acq_lat(i));
        pad      = zeros(1,5);
        pad(i+1) = floor(diff/2);
        coils    = utils.pad(coils, pad, 0, 'both');
        if nargout >= 3 && i ~= 1
            mask    = utils.pad(mask, pad([3 4]), 0, 'both');
        end
        if mod(diff,2)
            % odd crop
            pad(i+1) = 1;
            if mod(acq_lat(i),2)
                % odd lattice
                side = 'pre';
            else
                % even lattice
                side = 'post';
            end
            coils = utils.pad(coils, pad, 0, side);
            if nargout >= 3 && i ~= 1
                mask = utils.pad(mask, pad([3 4]), 0, side);
            end
        end
    end
end
% Update acquisition lattice
Nr = size(coils,2);
N1 = size(coils,3);
N2 = size(coils,4);
acq_lat = [Nr N1 N2];

% Inverse Fourier transform
coils = utils.ifft(coils, find(do_fft)+1);

% Crop dimensions where recon_fov < acq_fov
recon_lat = recon_lat ./ fov_factor;
for i=1:3
    if recon_lat(i) < acq_lat(i)
        crop   = ceil(acq_lat(i) - recon_lat(i));
        S      = struct;
        S.type = '()';
        if mod(acq_lat(i),2)
            % odd lattice
            if mod(crop,2)
                % odd crop
                idx = ceil(crop/2):(acq_lat(i)-ceil(crop/2));
            else
                % even crop
                idx = (crop/2+1):(acq_lat(i)-(crop/2));
            end
        else
            % even lattice
            if mod(crop,2)
                % odd crop
                idx = ceil(crop/2+1):(acq_lat(i)-floor(crop/2));
            else
                % even crop
                idx = (crop/2+1):(acq_lat(i)-(crop/2));
            end
        end
        S.subs = repmat({':'}, [1 5]);
        S.subs{i+1} = idx;
        coils  = subsref(coils, S);
        if nargout >= 3 && i ~= 1
            S.subs = S.subs([3 4]);
            mask = subsref(mask, S);
        end
    end
end

% Permute  [ch rd k1 k2 ct] -> [k1 k2 rd ch ct]
coils = permute(coils, [3 4 2 1 5]);
