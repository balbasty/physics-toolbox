function varargout = push(f, varargin)
% Push an image with a spatial transformation.
%
%   This function is a generic wraper of spm_diffeo('push') that can 
%   take deformation fields OR affine matrices as input.
%   Pushing is the adjoint operation of pulling. It shares the value of a
%   pushed voxels with all originally contributing locations.
%   Push assumes trilinear interpolation.
%
% FORMAT [imout,count] = utils.push(im, def, [dim])
%        [imout,count] = utils.push(im, mat, [dim])
% im    - Input image
% def   - Deformation field (voxel-to-voxel mapping)
% mat   - Affine transformation matrix (voxel-to-voxel mapping)
% dim   - Dimensions of the output image [default: same as input]
% imout - Deformed image
% count - Count (or contribution) image
%
% FORMAT utils.push(..., 'e')
%        utils.push(..., 'extrapolate')
% Extrapolate data outside of the original field of view [default: false]
%
% FORMAT utils.push(..., 'e')
%        utils.push(..., 'extrapolate')
% Extrapolate data outside of the original field of view [default: false]
%
% FORMAT utils.push(..., 'c')
%        utils.push(..., 'circulant')
%        utils.push(..., 'n')
%        utils.push(..., 'neumann')
% Use circulant ('c') or Neumann ('n') boundary conditions [default: 'n']

    % - Parse arguments
    y   = [];
    M   = [];
    dm  = {};
    ex  = false;
    bnd = 1;
    while ~isempty(varargin)
        if isequal(size(varargin{1}), [4 4])
        % Affine matrix
            M = varargin{1};
            if numel(varargin) > 1 && isnumeric(varargin{2})
                varargin = varargin(2:end);
                dm{1} = varargin{1};
            end
        elseif isnumeric(varargin{1})
        % Deformation field
            y = varargin{1};
            if numel(varargin) > 1 && isnumeric(varargin{2})
                varargin = varargin(2:end);
                dm{1} = varargin{1};
            end
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
        fun = 'pushc';
    else
        fun = 'push';
    end

    % Build transformation field
    if isempty(y)
        dmi = [size(f) 1];
        dmi = dmi(1:3);
        y = zeros(dmi, 'single');
        [y(:,:,:,1),y(:,:,:,2),y(:,:,:,3)] = ndgrid(1:dmi(1), 1:dmi(2), 1:dmi(3));
        y = reshape(y, [], 3);
        y = bsxfun(@plus, y * M(1:3,1:3)', M(1:3,4)');
        y = reshape(y, [dmi 3]);
    end

    % - Push
    spm_diffeo('boundary', bnd);
    [varargout{1:nargout}] = spm_diffeo(fun, single(f()), y, dm{:});

end