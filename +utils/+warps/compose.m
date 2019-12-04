function t = compose(varargin)
% Compose spatial transformations (warps or affine matrices)
%
% FORMAT t = compose(t1,t2,...,[dm])
% t<#> - Warping field OR affine matrix
% dm   - Output dimensions
% t    - Composed field or matrix: t = t1 o t2 o t3 ...

% -------------------------------------------------------------------------
% Check if output dimensions provided
dm = [];
if ~isempty(varargin) ...
        && size(varargin{end},1) == 1 ...
        && size(varargin{end},2) <= 3 ...
        && numel(size(varargin{end})) == 2
    dm = varargin{end};
    varargin = varargin(1:end-1);
end

% -------------------------------------------------------------------------
% Pass 1: collapse all affine matrices
j = 0;
for i=1:numel(varargin)
    if ismatrix(varargin{i})
        if j == 0
            j = i;
        else
            varargin{j} = compose_mat_mat(varargin{j}, varargin{i});
            varargin{i} = [];
        end
    else
        j = 0;
    end
end
varargin = varargin(~cellfun(@isempty, varargin));

% -------------------------------------------------------------------------
% Pass 2: perform all 'warp o affine' operations
j = 0;
for i=1:numel(varargin)
    if ~ismatrix(varargin{i})
        j = i;
    elseif j > 0
        varargin{j} = compose_warp_mat(varargin{j}, varargin{i});
        varargin{i} = [];
    end
end
varargin = varargin(~cellfun(@isempty, varargin));

% -------------------------------------------------------------------------
% Pass 3: create left-most warp if needed
if ismatrix(varargin{1})
    if isempty(dm)
        for j=2:numel(varargin)
            if ~ismatrix(varargin{j})
                dm = [size(varargin{j}) 1 1];
                dm = dm(1:3);
                brek
            end
        end
    end
    if ~isempty(dm)
        varargin{1} = compose_warp_mat(utils.warp.identity(dm), varargin{1});
    end
end

% -------------------------------------------------------------------------
% Pass 4: compose all warps
t = varargin{end};
varargin = varargin(1:end-1);
while ~isempty(varargin)
    t = compose_warp_warp(varargin{end}, t);
    varargin = varargin(1:end-1);
end

% =========================================================================
%   HELPERS
% =========================================================================

function M = compose_mat_mat(M1,M2)
M = M1 * M2;

function t = compose_warp_mat(t,M)
dm = [size(t) 1];
dm = dm(1:3);
t  = reshape(t, [], 3);
t  = bsxfun(@plus, t * M(1:3,1:3)', M(1:3,4)');
t  = reshape(t, [dm 3]);

function t = compose_warp_warp(t,t2)
spm_diffeo('boundary', 0);
dm = [size(t) 1];
dm = dm(1:3);
id = utils.warps.identity(dm);
t  = t - id;
t  = spm_diffeo('bsplins', single(t), single(t2), [1 1 1   1 1 1]);
t  = t + id;

