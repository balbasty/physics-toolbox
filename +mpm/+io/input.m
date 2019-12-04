function in = input(varargin)
% Create structure of input data.
%
%   This function extracts all necessary metadata from the input files and
%   organises them so that they are easy to process. The generated
%   structure can be used as an input to our fitting functions such as 
%   mpm.nonlin.fit, mpm.estatics.nonlin.fit or mpm.estatics.loglin.fit.
%   All input files are memory mapped to save RAM, unless the 'load' 
%   option is specified. Files should be nifti, and important metadata
%   (flip angle, TR, TE, ...) should be embedded in the 'description' field
%   of their header, or in a JSON sidecar file.
%
% FORMAT in = mpm.io.input(vol1, vol2, ..., ['load'|'map'])
% ------
%
% Each `vol` input should be, either:
%   . a (4D) nifti object or filename
%   . a cell of (3D) nifti objects or filenames
%   . eventually, each file or object can be accompanied with metadata 
%     embedded in a structure.
% 
% . if the 'load' option is specified, all volumes are loaded in memory.
% . if the 'map' option is specified, volumes are memory-mapped [default].
%
% EXAMPLES
% --------
% >> % One 4D volume with embedded metadata
% >> vol  = 'path/to/4Dvol.nii';
% >>
% >> % Several 3D volumes with embedded metadata
% >> vol  = {'path/to/echo1.nii'
%            'path/to/echo2.nii'
%            'path/to/echo3.nii'};
% >>
% >> % One 4D volume with specified metadata
% >> info = struct('FA', 6, ...
%                  'TR', 2.3, ...
%                  'echoes', struct('TE', {1.2, 2.4, 3.6}));
% >> vol  = {'path/to/4Dvol.nii', info};
% >>
% >> % Several 3D volumes with specified metadata
% >> vol  = {{'path/to/echo1.nii', struct('TE', 1.2)}
%            {'path/to/echo2.nii', struct('TE', 2.4)}
%            {'path/to/echo3.nii', struct('TE', 3.6)}
%            struct('FA', 6, 'TR', 2.3)};
%
% OUTPUT
% ------
% in{v}
%   .type       (PDw,T1w,MTw,B1+,B1-)
%   .seq        (SPGR,BS,...)
%   .dim        (volume dimensions)
%   .mat0       (original voxel-to-world matrix)
%   .mat        (registered voxel-to-world matrix)
%   .trf        (transformation matrix)
%   .nii(f)     (nifti file[s]) 
%   [if SPGR]
%   .field      (field strength in Tesla)
%   .FA         (nominal flip angle in radian)
%   .TR         (repetition time in second)
%   .SO         (magnetisation-transfer pulse)
%   .var        (noise variance: geometric mean across echoes)
%   .echoes{e}
%       .TE     (echo time in second)
%       .mean   (mean tissue intensity)
%       .var    (noise variance)
%       .dat    (file_array OR array [if 'load'])

% -------------------------------------------------------------------------
% Check load/map option
% -------------------------------------------------------------------------
load = false;
if ~isempty(varargin) ...
        && ischar(varargin{end}) ...
        && any(strcmpi(varargin{end}, {'load' 'map'}))
    load = strcmpi(varargin{end}, 'load');
    varargin = varargin(1:(end-1));
end

% -------------------------------------------------------------------------
% Create structure
% -------------------------------------------------------------------------

V  = numel(varargin);
in = cell(1,V);

for v=1:V
    fprintf('%d', v);
    arg1    = varargin{v};
    vol     = struct;
    % ---------------------------------------------------------------------
    % Fill everything that is provided (filenames, provided info)
    
    % -----------
    % First level
    if ischar(arg1) || isa(arg1, 'nifti')
        vol.nii = nifti(arg1);
    elseif iscell(arg1)
        vol.nii = nifti;
        % ------------
        % Second level
        if isstruct(arg1{end})
            volinfo = arg1{end};
            arg1 = arg1(1:(end-1));
            fields = fieldnames(volinfo);
            for f=1:numel(fields)
                vol.(fields{f}) = volinfo.(fields{f});
            end
        end
        N = numel(arg1);
        subinfo = cell(1,N);
        for n=1:N
            arg2 = arg1{n};
            if ischar(arg2) || isa(arg2, 'nifti')
                vol.nii(n) = nifti(arg2);
            elseif iscell(arg2)
                % -----------
                % Third level
                if isstruct(arg2{end})
                    subinfo{n} = arg2{end};
                    arg2 = arg2(1:(end-1));
                end
                arg3 = arg2{1};
                if ischar(arg3) || isa(arg3, 'nifti')
                    vol.nii(n) = nifti(arg3);
                else
                    error('Unsupported input type');
                end
                if numel(arg2) > 1
                    warning('Unsupported input type. Skipping it...');
                end
            else
                warning('Unsupported input type. Skipping it...');
            end
        end
    else
        error('Unsupported input type');
    end
    
    % ---------------------------------------------------------------------
    % Extract metadata
    if ~isfield(vol, 'dim')
        vol.dim  = [vol.nii(1).dat.dim 1];
        vol.dim  = vol.dim(1:3);
    end
    if ~isfield(vol, 'mat0')
        if strcmpi(vol.nii(1).mat0_intent, 'Scanner')
            vol.mat0 = vol.nii(1).mat0;
        elseif strcmpi(vol.nii(1).mat_intent, 'Scanner')
            vol.mat0 = vol.nii(1).mat;
        end
    end
    if ~isfield(vol, 'trf')
        vol.trf  = eye(4);
    end
    if ~isfield(vol, 'mat')
        vol.mat  = vol.trf\vol.mat0;
    end
    if ~isfield(vol, 'seq')
        vol.seq  = guessfield('seq', vol.nii);
    end
    if ~isfield(vol, 'type')
        vol.type = guessfield('type', vol.nii);
    end
    if strcmpi(vol.seq, 'SGE')
        if ~isfield(vol, 'field')
            vol.field = guessfield('field', vol.nii);
        end
        if ~isfield(vol, 'FA')
            vol.FA = guessfield('FA', vol.nii);
        end
        if ~isfield(vol, 'TR')
            vol.TR = guessfield('TR', vol.nii);
        end
        if ~isfield(vol, 'SO')
            vol.SO = guessfield('SO', vol.nii);
        end
        if numel(vol.nii) == 1
            N = size(vol.nii.dat, 4);
        else
            N = numel(vol.nii);
        end
        vol.echoes = cell(1,N);
        [vol.echoes{:}] = deal(struct);
        smean = 0;
        for n=1:N
            fprintf('.');
            if ~isfield(vol.echoes{n}, 'TE')
                vol.echoes{n}.TE = guessfield('TE', vol.nii, n);
            end
            % --- Setup data array
            vol.echoes{n}.dat = setupdat(vol.nii, n, load);
            if ~isfield(vol.echoes{n}, 'mean') || ~isfield(vol.echoes{n}, 'var')
                [s,mu] = utils.noise_estimate(vol.echoes{n}.dat());
                if ~isfield(vol.echoes{n}, 'mean')
                    vol.echoes{n}.mean = mu;
                end
                if ~isfield(vol.echoes{n}, 'var')
                    vol.echoes{n}.var = s;
                end
            end
            smean = smean + log(vol.echoes{n}.var);
        end
        smean = exp(smean/N);
        if ~isfield(vol, 'var')
            vol.var = smean;
        end
    end
    
    % ---------------------------------------------------------------------
    % Save current volume
    in{v} = vol;
    fprintf(' ');
end
fprintf('\n');

% -------------------------------------------------------------------------

function dat = setupdat(nii, n, load)
% Generate a 3D array (or file_array) from a (3D or 4D) nifti object.
%
% FORMAT dat = setupdat(nii, n, load)
% nii  - Nifti structure[s]
% n    - Index of array to setup
% load - If true, load in memory (else, memory map using file_array)
% dat  - 3D array or (read only) file_array

if numel(nii) == 1
    % Create a 3D file_array out of a 4D file_array
    dat = nii.dat;
    dtype = lower(dat.dtype);
    dtype = strsplit(dtype, '-');
    dtype = dtype{1};
    switch dtype
        case {'int8' 'uint8'}
            nbytes = 1;
        case {'int16' 'uint16'}
            nbytes = 2;
        case {'rgb24'}
            nbytes = 3;
        case {'int32' 'uint32' 'float32'}
            nbytes = 4;
        case {'int64' 'uint64' 'float64' 'complex64'}
            nbytes = 8;
        case {'float128' 'complex128'}
            nbytes = 16;
        case {'complex256'}
            nbytes = 32;
    end
    dim = [dat.size 1];
    dim = dim(1:3);
    dat.offset = dat.offset + prod(dim) * nbytes * (n-1);
    dat.dim = dim;
else
    dat = nii(n).dat;
end
dat.permission = 'ro';
if load
    dat = dat();
end

% -------------------------------------------------------------------------

function val = guessfield(field, nii, n)
% Try to guess parameter value from filename/nifti header/json sidecar.
%
% FORMAT val = guessfield(field, nii, n)
% field - Parameter name
% nii   - Nifti object[s]
% n     - Index of 3D volume [default: 1]
% val   - Value [NaN if not found]

if nargin < 3
    n = 1;
end
if numel(nii) > 1
    nii = nii(n);
end
fname   = nii.dat.fname;
descrip = nii.descrip;

val = from_fname(field, fname); % seq/type
if ~isnan(val)
    return
end
val = from_descrip(field, descrip); % TE/TR/FA/SO/field
if ~isnan(val)
    return
end
val = from_json(field, fname, n); % seq/type/TE/TR/FA/SO/field

% -------------------------------------------------------------------------

function prm = from_fname(type, fname)

prm = NaN;
[path,fname] = fileparts(fname);
path = strsplit(path, filesep);

switch lower(type)
    case 'type'
        if numel(fname) >= 3
            switch lower(fname(1:3))
                case 't1w'
                    prm = 'T1w';
                case 'pdw'
                    prm = 'PDw';
                case 'mtw'
                    prm = 'MTw';
            end
            if ~isnan(prm), return; end
        end
        if ~isempty(path) && numel(path{end}) >= 3
            switch lower(path{end}(1:3))
                case 't1w'
                    prm = 'T1w';
                case 'pdw'
                    prm = 'PDw';
                case 'mtw'
                    prm = 'MTw';
            end
            if ~isnan(prm), return; end
        end
    case 'seq'
        % spoiled gradient echo (SGE)
        if ~isempty(regexpi(fname, 'flash'))
            prm = 'SGE';
            return
        end
        if ~isempty(path) && ~isempty(regexpi(path{end}, 'flash'))
            prm = 'SGE';
            return
        end
        % spin echo stimulated echo (SESTE)
        if ~isempty(regexpi(fname, 'seste'))
            prm = 'SESTE';
            return
        end
        if ~isempty(path) && ~isempty(regexpi(path{end}, 'seste'))
            prm = 'SESTE';
            return
        end
end

% -------------------------------------------------------------------------

function prm = from_descrip(type, descrip)

prm = NaN;

switch lower(type)
    case {'te' 'tr' 'fa' 'so'}
        res = regexp(descrip, [upper(type) '=(?<prm>([\d.]+|no|MT))(?<unit>([m]?s|deg)?)[ /]?'],'names');
        if ~isempty(res)
            switch lower(type)
                case {'te' 'tr' 'fa'}
                    prm = str2double(res.prm);
                case 'so'
                    prm = upper(res.prm);
            end
            switch res.unit
                case 'ms'
                    prm = prm * 1E-3;
                case 'deg'
                    prm = prm * pi / 180;
            end
        end
    case {'field'}
        res = regexp(descrip, '^(?<prm>([\d.]+))T ','names');
        if ~isempty(res)
            prm = str2double(res.prm);
        end
end

% -------------------------------------------------------------------------

function prm = from_json(type, fname, n)

if nargin < 3, n = 1; end

prm = NaN;

switch lower(type)
    case 'te'
        field = 'EchoTime';
    case 'tr'
        field = 'RepetitionTime';
    case 'fa'
        field = 'FlipAngle';
    case 'so'
        field = 'ScanOptions';
    case 'seq'
        field = 'SequenceName';
    otherwise
       return
end


[dir,name] = fileparts(fname);
fname = fullfile(dir, [name '.json']);

if exist(fname, 'file')
    meta = spm_jsonread(fname);
    if isfield(meta, 'acqpar') && isfield(meta.acqpar, field)
        prm = meta.acqpar.(field);
    end
    switch lower(type)
        case {'te' 'tr'}
            prm = str2double(res.prm);
            prm = prm * 1E-3;
        case 'fa'
            prm = str2double(res.prm);
            prm = prm * pi / 180;
        case 'so'
            prm = upper(prm);
            % prm = strcmpi(prm, 'mt');
        case 'seq'
            prm = ~isempty(regexpi(prm, '(fl3d|flash)'));
            if prm, prm = 'SGE';
            else,   prm = NaN; end
    end
end