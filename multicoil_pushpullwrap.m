function x = multicoil_pushpullwrap(x, msk, dirs)
% FORMAT y = multicoil_pushpullwrap(x,msk,dir)
% x   - Non-wrapped image [n1 n2 n3 n4 1/2 n6]
% msk - Mask of the k-space sampling scheme [n1 n2 (n3)]
% dir - List of accelerated directions (default: from msk)
% y   - Pulled image [n1 n2 n3 n4 1/2 n6]

if isempty(msk)
    return
end
if nargin < 3
    dirs = logical([0 0 0]);
    mskdim = [size(msk) 1];
    dirs(1) = mskdim(1) == 1 || all(all(all((msk(1:2:2*floor(mskdim(1)/2),:,:) - msk(2:2:2*floor(mskdim(1)/2),:,:))==0,1),2),3);
    dirs(2) = mskdim(2) == 1 || all(all(all((msk(:,1:2:2*floor(mskdim(2)/2),:) - msk(:,2:2:2*floor(mskdim(2)/2),:))==0,1),2),3);
    dirs(3) = mskdim(3) == 1 || all(all(all((msk(:,:,1:2:2*floor(mskdim(3)/2)) - msk(:,:,2:2:2*floor(mskdim(3)/2)))==0,1),2),3);
    dirs = find(~dirs);
end
dirs = dirs(:)';
if isempty(dirs)
    return
end

msk2d = numel(size(msk)) == 2;
msk   = reshape(msk, [], 1);

for m=1:size(x,6) % echoes/contrasts
for n=1:size(x,4) % coils/channels

    % Read image
    if size(x,5) == 2
        x1 = x(:,:,:,n,1,m) + 1i*x(:,:,:,n,2,m);
    else
        x1 = x(:,:,:,n,1,m);
    end
    dim = [size(x1) 1 1];
    
    % FFT transform
    for dir=dirs
        x1 = ifftshift(x1,dir);
        x1 = fft(x1,[],dir);
        x1 = fftshift(x1,dir);
    end

    % Decimate k-space
    x1 = reshape(x1, numel(msk), []);
    if msk2d
        x1(~msk,:) = 0;
    else
        x1(~msk) = 0;
    end
    x1 = reshape(x1, dim);

    % Inverse FFT
    for dir=dirs
        x1 = ifftshift(x1,dir);
        x1 = ifft(x1,[],dir);
        x1 = fftshift(x1,dir);
    end

    % Save image
    if size(x,5) == 2
        x1 = cat(5,real(x1),imag(x1));
    end
    x(:,:,:,n,:,m) = x1;
    
end
end