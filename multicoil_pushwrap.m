function x = multicoil_pushwrap(x, msk, dirs, dom)
% FORMAT y = multicoil_pushwrap(x,smp,[dir],[dom])
% x   - Wrapped image [n1 n2 n3 n4 n5]
% smp - Mask of the k-space sampling scheme [n1 n2 (n3)]
% dir - List of accelerated directions (default: from msk)
% dom - Input and output domains can be image ('im') or frequency ('freq')
%       [{'im' 'im'}]
% y   - Pushed image [n1 n2 n3 n4 n6]
%
% This operator unwraps an image according to a k-space sampling scheme. 
% This does not perform an actual unwrapping (it does not create a nice 
% object image) but merely performs the adjoint of the "wraping"
% operations, which we will name "pushing".
%
% Note:
% . n1/n2 = accelerated directions
% . n3    = readout direction
% . n4    = number of coils/channels
% . n5    = number of echoes/contrasts

if nargin < 4
    dom = {'im' 'im'};
end
if isempty(msk) && strcmpi(dom{1}, dom{2})
    return
end
if nargin < 3 || isempty(dirs)
    % trying to detect accelerated directions automatically based on the
    % shape of the sampling mask.
    dirs = logical([0 0 0]);
    mskdim = [size(msk) 1];
    dirs(1) = mskdim(1) == 1 || all(all(all((msk(1:2:2*floor(mskdim(1)/2),:,:) - msk(2:2:2*floor(mskdim(1)/2),:,:))==0,1),2),3);
    dirs(2) = mskdim(2) == 1 || all(all(all((msk(:,1:2:2*floor(mskdim(2)/2),:) - msk(:,2:2:2*floor(mskdim(2)/2),:))==0,1),2),3);
    dirs(3) = mskdim(3) == 1 || all(all(all((msk(:,:,1:2:2*floor(mskdim(3)/2)) - msk(:,:,2:2:2*floor(mskdim(3)/2)))==0,1),2),3);
    dirs = find(~dirs);
end
dirs = dirs(:)';
if isempty(dirs) && strcmpi(dom{1}, dom{2})
    return
end

gpu_on = isa(msk, 'gpuArray');
if gpu_on, loadarray = @loadarray_gpu;
else,      loadarray = @loadarray_cpu; end

msk    = reshape(msk, [], 1);
mskdim = numel(size(msk));

for m=1:size(x,5) % echoes/contrasts
for n=1:size(x,4) % coils/channels

    % Read image
    x1 = loadarray(x(:,:,:,n,m), @single);
    dim = [size(x1) 1];
    dim = dim(1:3);
    
    % Allocate output imeag
    y1 = zeros([numel(msk) dim(mskdim+1:end)], 'single');
    if gpu_on
        y1 = gpuArray(y1);
    end
    
    if strcmpi(dom{1}, 'im')
        % Inverse FFT
        for dir=dirs
            x1 = ifftshift(x1,dir);
            x1 = ifft(x1,[],dir);
            x1 = fftshift(x1,dir);
        end
        
        % Subsample k-space
        x1 = reshape(x1, numel(msk), []);
        x1 = x1(msk,:);
    end

    y1(msk,:) = x1;

    % FFT
    if strcmpi(dom{2}, 'im')
        for dir=dirs
            x1 = ifftshift(x1,dir);
            x1 = fft(x1,[],dir);
            x1 = fftshift(x1,dir);
        end
    end

    % Save image
    x(:,:,:,n,m) = x1;
    
end
end