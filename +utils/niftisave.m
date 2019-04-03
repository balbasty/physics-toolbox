function niftisave(fname, vol, varargin)
% Save a matlab array as Nifti file
%
% FORMAT niftisave(fname, vol, ...)
%
% REQUIRED
% --------
% fname   - File name
% vol     - Matlab array
%
% KEYWORD
% -------
% mat     - Voxel to world matrix   [diag(vs)]
% vs      - Voxel size              [from mat or [1 1 1]]
% dtype   - file_array data type    [from array type]
% slope   - nifti data slope        [1]
% inter   - nifti data intercept    [0]
% descrip - Description             ['']
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Parse arguments
% -------------------------------------------------------------------------

p = inputParser;
p.FunctionName = 'niftisave';
p.addRequired('fname',            @ischar);
p.addRequired('vol',              @isnumeric);
p.addParameter('mat',     eye(4), @ismat);
p.addParameter('vs',      [],     @isnumeric);
p.addParameter('dtype',   '',     @istype);
p.addParameter('slope',   1,      @isnumericscalar);
p.addParameter('inter',   0,      @isnumericscalar);
p.addParameter('descrip', '',     @ischar);
p.parse(fname, vol, varargin{:});
mat     = p.Results.mat;
vs      = p.Results.vs;
dtype   = p.Results.dtype;
slope   = p.Results.slope;
inter   = p.Results.inter;
descrip = p.Results.descrip;


% -------------------------------------------------------------------------
% Process arguments
% -------------------------------------------------------------------------

vol = gather(vol);

if isempty(dtype)
    dtype = map_type(class(vol),isreal(vol));
end

if ~isempty(vs)
    vs0 = sqrt(sum(mat(1:3,1:3).^2));
    mat(1:3,1:3) = mat(1:3,1:3) * diag(vs./vs0);
end

% -------------------------------------------------------------------------
% Create nifti file
% -------------------------------------------------------------------------

dat         = file_array(fname, size(vol), dtype, 0, slope, inter);
nii         = nifti;
nii.dat     = dat;
nii.mat     = mat;
nii.mat0    = mat;
nii.descrip = descrip;
create(nii);
nii.dat(:)  = vol(:);

end

% =========================================================================
% Helper functions
% =========================================================================

function ok = ismat(x)
    ok = isnumeric(x) && issame(size(x), [4 4]);
end

function ok = istype(x)
    list_types = {'binary' ...
        'float32' 'float64' 'float128' ...
        'complex32' 'complex64' 'complex128' ...
        'uint8' 'uint16' 'uint32' 'uint64' ...
        'int8' 'int16' 'int32' 'int64'};
    ok = ischar(x);
    x = strsplit(x, '-');
    ok = ok && (isempty(x{1}) || any(strcmpi(x{1},list_types)));
    if numel(x) > 1
        ok = ok && any(strcmpi(x{2}, {'le' 'be'}));
    end
end

function ok = isnumericscalar(x)
    ok = isnumeric(x) && isscalar(x);
end

function dtype = map_type(mtype,isreal)
    if isreal
        switch lower(mtype)
            case 'single'
                dtype = 'float32';
            case 'double'
                dtype = 'float64';
            case 'uint8'
                dtype = 'uint8';
            case 'uint16'
                dtype = 'uint16';
            case 'uint32'
                dtype = 'uint32';
            case 'uint64'
                dtype = 'uint64';
            case 'int8'
                dtype = 'int8';
            case 'int16'
                dtype = 'int16';
            case 'int32'
                dtype = 'int32';
            case 'int64'
                dtype = 'int64';
            case 'logical'
                dtype = 'binary';
            otherwise
                warning('Real data type %s not handled. Using float32 instead.', mtype);
                dtype = 'float32';
        end
    else
        switch lower(mtype)
            case 'single'
                dtype = 'complex64';
            case 'double'
                dtype = 'complex128';
            otherwise
                warning('Complex data type %s not handled. Using complex64 instead.', mtype);
                dtype = 'complex64';
        end
    end
end