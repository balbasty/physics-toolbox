function o = opt(o)
% Set default options.
%
% FORMAT opt_out = mpm.estatics.nonlin.opt(opt_in)
% opt_in  - User-defined option structure [empty]
% opt_out - All user-defined options are kept untouched, the other are set 
%           to their default value.

% --- Default values
if nargin < 1, o = struct; end
o = setdefault(o, 'nbscales',         5);         % Number of scales
o = setdefault(o, 'nbiter',           3);         % Number of iterations per scale
o = setdefault(o, 'tolerance',        1E-4);      % Gain threshold for early stopping (per scale)
o = setdefault(o, 'out.folder',       '.');       % Output folder
o = setdefault(o, 'out.fname',        '.nii');    % Suffix for output files 
o = setdefault(o, 'out.mem',          'map');     % map/load output volumes
o = setdefault(o, 'reg.mode.default', [0 1]);     % Absolute/Membrane regul. (0=None|1=L1|2=L2)
o = setdefault(o, 'reg.prec.default', [1E0 5E3]); % Absolute/Membrane precision
o = setdefault(o, 'reg.prec.R2s',     [1E0 5E0]);
o = setdefault(o, 'reg.mean.default', NaN);       % Mean (NaN = from minilogfit)
o = setdefault(o, 'reg.uncertainty' , 1E-3);      % RLS smoother (value|'bayes')
o = setdefault(o, 'vs',               NaN);       % Reconstruction voxel size (Nan=from input)
o = setdefault(o, 'fov',              0);         % Field of view (0=bounding box|n=index of input volume)
o = setdefault(o, 'coreg',            true);      % Co-register volumes first
o = setdefault(o, 'init',             'mean');    % Initialisation mode: ('mean'|'logfit'|'minilogfit')
o = setdefault(o, 'subsample',        Inf);       % Subsampling distance (Inf=no subsampling)
o = setdefault(o, 'threads',          -1);        % Number of threads (-1=all)
o = setdefault(o, 'verbose',          1);         % Verbosity (0=quiet|[1]=print|2=plot)
o = setdefault(o, 'solver.type',      'relax');   % Solver type ('relax'|'cg')
o = setdefault(o, 'solver.nbiter',    10);        % Number of iterations of the linear solver
o = setdefault(o, 'solver.tolerance', 0);         % Solver gain threshold 
o = setdefault(o, 'solver.verbose',   true);      % Solver verbosity 
o = setdefault(o, 'solver.sumtype',   'double');  % Solver accumulator type ('double'|'native')
o = setdefault(o, 'solver.precond',   true);      % [CG only] use preconditioner 

% --- Reformat options
o.vs = utils.pad(o.vs(:)', [0 3-numel(o.vs)], 'replicate', 'post');

% --- Get solver handle
switch lower(o.solver.type)
    case 'relax'
        o.solver.fun = @mpm.l1.solve.relax;
    case 'cg'
        o.solver.fun = @mpm.l1.solve.cg;
end

% --- Ensure output folder exists
path = strsplit(o.out.folder, filesep);
if isempty(path{1}), path{1} = filesep; end
% TODO: ^ not sure about this, especially on windows...
for i=1:numel(path)
    folder = fullfile(path{1:i});
    if ~exist(folder, 'dir')
        st = mkdir(folder);
        if ~st
            error('Cannot create output folder %s', folder);
        end
    end
end

% -------------------------------------------------------------------------
function o = setdefault(o, field, value, exists)
% Set default value in option structure. The value is only set if none 
% existed previously.
%
% FORMAT opt = setdefault(opt, field, value)
% opt   - Structure
% field - Hierarchy if fields (cell, or '.' separated)
% value - Default value
%
% EXAMPLE
% >> opt.my.field1 = 1;
% >> opt = setdefault(opt, 'my.field1', 2);
% >> opt.my.field1
% ans = 1
% >> opt = setdefault(opt, 'my.field2', 3);
% >> opt.my.field2
% ans = 3
if nargin < 4, exists = false; end
if ~iscell(field), field = strsplit(field, '.'); end

if isempty(field)
    if ~exists, o = value; end
else
    exists = isfield(o, field{1});
    if ~exists, o.(field{1}) = []; end
    o.(field{1}) = setdefault(o.(field{1}), field(2:end), value, exists);
end