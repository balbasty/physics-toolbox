function x = multicoil_pullwrap(x, msk, dirs)
% FORMAT y = multicoil_pullwrap(x,msk,[dir],[dom])
% x   - Non-wrapped image [n1 n2 n3 n4 n6]
% msk - Mask of the k-space sampling scheme [n1 n2 (n3)]
% dir - List of accelerated directions (default: from msk)
% dom - Input and output domains can be image ('im') or frequency ('freq')
%       [{'im' 'im'}]
% y   - Pulled image [n1 n2 n3 n4 n6]
%
% This operator wraps an image according to a k-space sampling scheme.
%
% Note:
% . n1/n2 = accelerated directions
% . n3    = readout direction
% . n4    = number of coils/channels
% . n5    = number of echoes/contrasts

if isempty(msk) && strcmpi(dom{1}, dom{2})
    return
end
if nargin < 4
    dom = {'im' 'im'};
end
if nargin < 3 || isempty(dirs)
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
    
    if strcmpi(dom{1}, 'im')
        % FFT
        for dir=dirs
            x1 = ifftshift(x1,dir);
            x1 = fft(x1,[],dir);
            x1 = fftshift(x1,dir);
        end
    end

    if strcmpi(dom{2}, 'im')
        % Decimate k-space
        x1 = reshape(x1, numel(msk), []);
        x1(~msk,:) = 0;
        x1 = reshape(x1, dim);
        
        % Inverse FFT
        for dir=dirs
            x1 = ifftshift(x1,dir);
            x1 = fft(x1,[],dir);
            x1 = fftshift(x1,dir);
        end
    else
        % Undersample k-space
        x1 = reshape(x1, numel(msk), []);
        x1 = x1(msk,:);
        x1 = reshape(x1, [size(x1,1) dim(mskdim+1:end)]);
    end

    % Save image
    x(:,:,:,n,m) = x1;
    
end
end