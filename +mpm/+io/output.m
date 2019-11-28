function out = output(prefix, dim, mat, val, typ, opt)
% Create structure for holding output (model) data.
%
% FORMAT out = mpm.io.output(prefix, dim, [mat], [val], [typ], [opt])
% prefix {N}     - Cell of unique names/types/prefixes
% dim    {3 N}   - Dimensions of each volume
% mat    {4 4 N} - Orientation matrix of each volume [eye(4)]
% val    {N}     - Initial value of each volume [0]
% typ    {N}     - Data type of each volume ['single']
% opt            - Structure of optional parameters with fields:
%   . folder - Output folder
%   . fname  - Suffix (with extension) of each file
%   . mem    - 'map'/'load'
%
% Note that, if opt.mem == 'load', the nifti file is not created.
% Use `mpm.io.create(out)` to create the nifti file and write the volumes.

out  = struct;

if nargin < 6, opt = struct; end
if ~isfield(opt, 'folder'), opt.folder = '.';    end
if ~isfield(opt, 'fname'),  opt.fname  = '.nii'; end
if ~isfield(opt, 'mem'),    opt.mem    = 'map';  end
if nargin < 5 || isempty(typ), typ = {'single'}; end
if nargin < 4 || isempty(val), val = 0;          end
if nargin < 3 || isempty(mat), mat = eye(4);     end

N = numel(prefix);
if size(dim,1) == 1, dim = dim(:); end
dim = utils.pad(dim, [0 N-size(dim,2)], 'replicate', 'post');
mat = utils.pad(mat, [0 0 N-size(mat,3)], 'replicate', 'post');
val = utils.pad(val(:), [N-numel(val) 0], 'replicate', 'post');
if ~iscell(typ), typ = {typ}; end
typ = utils.pad(typ(:), [N-numel(typ) 0], 'replicate', 'post');

if opt.fname(1) ~= '.' && opt.fname(1) ~= '_'
    opt.fname = ['_' opt.fname];
end

% -------------------------------------------------------------------------
% Loop over output volumes
% -------------------------------------------------------------------------
for i=1:numel(prefix)
    pre   = prefix{i};
    fname = fullfile(opt.folder, [prefix{i} opt.fname]);
    if strcmpi(opt.mem, 'map')
        out.(pre)        = nifti;
        out.(pre).dat    = file_array(fname, dim(:,i)', type_matlab_to_nifti(typ{i}));
        out.(pre).dat(:) = val(i);
        out.(pre).mat    = mat(:,:,i);
        create(out.(pre));
    else
        out.(pre).fname  = fname;
        out.(pre).dat    = zeros(dim(:,i)', typ{i}) + val(i);
        out.(pre).mat    = mat(:,:,i);
    end
end

% -------------------------------------------------------------------------
% Global fields
% -------------------------------------------------------------------------
if all(min(dim(1:3,:),[],2) == max(dim(1:3,:),[],2))
    out.dim = dim(1:3,1)';
    out.mat = mat(:,:,1);
end

function dtype = type_matlab_to_nifti(dtype,isreal)
if nargin < 2
    isreal = true;
end
switch lower(dtype)
    case 'single'
        if isreal
            dtype = 'float32';
        else
            dtype = 'complex64';
        end
    case 'double'
        if isreal
            dtype = 'float64';
        else
            dtype = 'complex128';
        end
end