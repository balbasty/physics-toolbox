%==========================================================================
function D = imgrad(X,vs,type) 
% Calculate spatial gradient of an image (with voxel size)
%
% FORMAT grad = utils.imgrad(img,(vs),(type))
% img   - Image (Nx * Nz * Nz * ...)
% vs    - voxel size [1]
% type  - Finite difference type '+'/'-'/['+-']/'-+'
% grad  - Gradients in all direction (Nx * Nz * Nz * ... * Ndim * Ntype)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% -------------------------------------------------------------------------
% Default parameters
if nargin < 2, vs = ones(1,'like',X); end
if nargin < 3, type = '+-'; end

dim    = size(X);
if numel(dim) == 2
    dim(3) = 1;
end
dimout = [prod(dim) numel(dim) numel(type) ];
vs     = utils.pad(vs(:)', [0 max(0,3-numel(dim))], 'replicate', 'post');

% -------------------------------------------------------------------------
% Compute each derivative that is asked for.
% The results of diff is of length N-1, so we pad the result with a zero.
D = zeros(dimout, 'like', X);
Z = sqrt(numel(type));
for i=1:numel(dim)
    if dim(i)  > 1
        dimpad    = zeros(size(dim));
        dimpad(i) = 1;
        diffX     = diff(X,1,i) ./ (vs(i) * Z);
        for t=1:numel(type)
            switch type(t)
                case '+'
                    D(:,i,t) = reshape(utils.pad(diffX, dimpad, 0, 'post'), [], 1);
                case '-'
                    D(:,i,t) = reshape(utils.pad(diffX, dimpad, 0, 'pre'), [], 1);
            end
        end
    end
end

% -------------------------------------------------------------------------
% Reorder dimensions > Nx * Nz * Nz * ... * Ndim * Ntype
% D = permute(D, [3 2 1]);
D = reshape(D, [dim numel(dim) numel(type)]);
end
%==========================================================================

% function G = imgrad(im,vx) 
% % Compute image gradient
% %__________________________________________________________________________
% % Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   
% 
% if nargin < 2, vx = ones([1 3],'like',im); end
% 
% Gx = cat(2, diff(im,1,2), zeros(size(im,1),1,size(im,3), 'like', im))./vx(1);
% Gy = cat(1, diff(im,1,1), zeros(1,size(im,2),size(im,3), 'like', im))./vx(2);
% Gz = cat(3, diff(im,1,3), zeros(size(im,1),size(im,2),1, 'like', im))./vx(3);  
% 
% G = cat(4,Gx,Gy,Gz);
% %==========================================================================