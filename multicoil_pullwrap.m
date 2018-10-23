function x = multicoil_pullwrap(x, msk, dirs)
% FORMAT y = multicoil_pullwrap(x,msk, dir)
% x   - Non-wrapped image [n1 n2 n3 n4]
% msk - Mask of the k-space sampling scheme [n1 n2 n3]
% dir - List of accelerated directions (default: [1 2])
% y   - Pulled image [n1 n2 n3 n4]
%
% This operator wraps an image according to a k-space sampling scheme

if nargin < 3
    dirs = [1 2];
end
dirs = dirs(:)';

msk = reshape(msk, [], 1);

for n=1:size(x,4)

    % Read image
    if size(x,5) == 2
        x1 = x(:,:,:,n,1) + 1i*x(:,:,:,n,2);
    else
        x1 = x(:,:,:,n);
    end
    dim = [size(x1) 1 1];
    
    % FFT transform
    for dir=dirs
        x1 = fft(x1,[],dir);
        x1 = fftshift(x1,dir);
    end

    % Decimate k-space
    x1 = reshape(x1, numel(msk), []);
    x1(~msk) = 0;
    x1 = reshape(x1, dim);

    % Inverse FFT
    for dir=dirs
        x1 = fftshift(x1,dir);
        x1 = ifft(x1,[],dir);
    end

    % Save image
    if size(x,5) == 2
        x1 = cat(5,real(x1),imag(x1));
    end
    x(:,:,:,n,:) = x1;
    
end