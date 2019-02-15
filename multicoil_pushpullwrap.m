function x = multicoil_pushpullwrap(x, msk, dirs, dom)
% FORMAT y = multicoil_pushpullwrap(x,msk,[dir],[dom])
% x   - Non-wrapped image [n1 n2 n3 n4 n5]
% msk - Mask of the k-space sampling scheme [n1 n2 (n3)]
% dir - List of accelerated directions (default: from msk)
% dom - Input and output domains can be image ('im') or frequency ('freq')
%       [{'im' 'im'}]
% y   - Pulled image [n1 n2 n3 n4 n5]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 4
    dom = {'im' 'im'};
end
if isempty(msk) && strcmpi(dom{1},dom{2})
    return
end
if nargin < 3 || isempty(dirs)
    dirs = logical([0 0 0]);
    mskdim = [size(msk) 1];
    dirs(1) = mskdim(1) == 1 || all(all(all(diff(msk,1)==0,1),2),3);
    dirs(2) = mskdim(2) == 1 || all(all(all(diff(msk,2)==0,1),2),3);
    dirs(3) = mskdim(3) == 1 || all(all(all(diff(msk,3)==0,1),2),3);
    dirs = find(~dirs);
end
dirs = dirs(:)';
if isempty(dirs) && strcmpi(dom{1},dom{2})
    return
end

gpu_on = isa(msk, 'gpuArray');
if gpu_on, loadarray = @loadarray_gpu;
else,      loadarray = @loadarray_cpu; end
gpu_in = isa(x, 'gpuArray');

invmsk = reshape(~msk, [], 1);

if isa(x, 'file_array')
    % Memory management seems to matters and load volumes one at a time
    for m=1:size(x,5) % echoes/contrasts
    for n=1:size(x,4) % coils/channels
        x(:,:,:,n,m) = pushpull(loadarray(x(:,:,:,n,m)), invmsk, dirs, dom);
    end
    end
else
    % The volume is already on memory, so we can process it at once
    x = pushpull(loadarray(x), invmsk, dirs, dom);
end
if ~gpu_in
    x = gather(x);
end


function x = pushpull(x, invmsk, dirs, dom)

% Performance note: On GPU, circshift is super slow. I did not find a nice
% workaround, so better not to do anything on GPU for the moment.
% Note that two circshift could be avoided by storing `ifftshift(mask)`
% instead of `mask`. The remaining two circshift could be removed if we
% used circulant boundary conditions rather than Neuman's for the
% sensitivity fields, in which case we could always work with shifted
% volumes.

dim = size(x);

if strcmpi(dom{1}, 'im')
    % FFT
    for dir=dirs
        x = ifftshift(x,dir);
        x = fft(x,[],dir);
        x = fftshift(x,dir);
    end
end

% Decimate k-space
x = reshape(x, numel(invmsk), []);
x(invmsk,:) = 0;
x = reshape(x, dim);

if strcmpi(dom{1}, 'im')
    % Inverse FFT
    for dir=dirs
        x = ifftshift(x,dir);
        x = ifft(x,[],dir);
        x = fftshift(x,dir);
    end
end
