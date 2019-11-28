function f = pull(f, varargin)
% Pull an image with a spatial transformation.
%
%   This function is a generic wraper of spm_diffeo('pull') that can 
%   take deformation fields OR affine matrices as input.
%   Pull assumes trilinear interpolation.
%
% FORMAT imout = utils.pull(im, def)
%        imout = utils.pull(im, mat, [dim])
% im    - Input image
% def   - Deformation field (voxel-to-voxel mapping)
% mat   - Affine transformation matrix (voxel-to-voxel mapping)
% dim   - Dimensions of the output image [default: same as input]
% imout - Deformed image
%
% FORMAT utils.pull(..., 'e')
%        utils.pull(..., 'extrapolate')
% Extrapolate data outside of the original field of view [default: false]
%
% FORMAT utils.pull(..., 'c')
%        utils.pull(..., 'circulant')
%        utils.pull(..., 'n')
%        utils.pull(..., 'neumann')
% Use circulant ('c') or Neumann ('n') boundary conditions [default: 'n']

    % - Parse arguments
    y   = [];
    M   = [];
    dm  = [];
    ex  = false;
    bnd = 1;
    while ~isempty(varargin)
        if isequal(size(varargin{1}), [4 4])
        % Affine matrix
            M = varargin{1};
            if numel(varargin) > 1 && isnumeric(varargin{2})
                varargin = varargin(2:end);
                dm = varargin{1};
            end
        elseif isnumeric(varargin{1})
        % Deformation field
            y = varargin{1};
        elseif ischar(varargin{1})
        % Option
            switch lower(varargin{1})
                case {'e' 'extrapolate'}
                    ex = true;
                case {'c' 'circulant'}
                    bnd = 0;
                case {'n' 'neumann'}
                    bnd = 1;
            end
        end
        varargin = varargin(2:end);
    end
    
    % - Select pulling function
    if ex
        fun = 'pullc';
    else
        fun = 'pull';
    end

    % - Affine case: build transformation field
    if isempty(y)
        if isempty(dm)
            dm = [size(f) 1];
            dm = dm(1:3);
        end
        y = zeros(dm, 'single');
        [y(:,:,:,1),y(:,:,:,2),y(:,:,:,3)] = ndgrid(1:dm(1), 1:dm(2), 1:dm(3));
        y = reshape(y, [], 3);
        y = bsxfun(@plus, y * M(1:3,1:3)', M(1:3,4)');
        y = reshape(y, [dm 3]);
    end

    % - Sample
    spm_diffeo('boundary', bnd);
    f = spm_diffeo(fun, single(f()), y);

end