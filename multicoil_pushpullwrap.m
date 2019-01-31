function x = multicoil_pushpullwrap(x, msk, dirs, dom)
% FORMAT y = multicoil_pushpullwrap(x,msk,[dir],[dom])
% x   - Non-wrapped image [n1 n2 n3 n4 n5]
% msk - Mask of the k-space sampling scheme [n1 n2 (n3)]
% dir - List of accelerated directions (default: from msk)
% dom - Input and output domains can be image ('im') or frequency ('freq')
%       [{'im' 'im'}]
% y   - Pulled image [n1 n2 n3 n4 n5]

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

msk = reshape(msk, [], 1);

for m=1:size(x,5) % echoes/contrasts
for n=1:size(x,4) % coils/channels

    % Read image
    x1 = loadarray(x(:,:,:,n,m));
    dim = [size(x1) 1 1];
    
    if strcmpi(dom{1}, 'im')
        % FFT
        for dir=dirs
            x1 = ifftshift(x1,dir);
            x1 = fft(x1,[],dir);
            x1 = fftshift(x1,dir);
        end
    end

    % Decimate k-space
    x1 = reshape(x1, numel(msk), []);
    x1(~msk,:) = 0;
    x1 = reshape(x1, dim);

    if strcmpi(dom{1}, 'im')
        % Inverse FFT
        for dir=dirs
            x1 = ifftshift(x1,dir);
            x1 = ifft(x1,[],dir);
            x1 = fftshift(x1,dir);
        end
    end

    % Save image
    x(:,:,:,n,m) = x1;
    
end
end