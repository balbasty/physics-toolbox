function x = multicoil_pushwrap(x, msk, dirs)
% FORMAT y = multicoil_pushwrap(x,msk)
% x   - Wrapped image [n1 n2 n3 n4]
% msk - Mask of the k-space sampling scheme [n1 n2 n3]
% dir - List of accelerated directions (default: [1 2])
% y   - Pushed image [n1 n2 n3 n4]
%
% This operator unwraps an image according to a k-space sampling scheme. 
% This does not perform an actual unwrapping (it does not create a nice 
% object image) but merely performs the transpose of the "wraping"
% operations, which we will name "pushing".

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
    
    % Inverse FFT
    % x1 = fftshift(fftshift(fftshift(x1,1),2),3);
    for dir=dirs
        x1 = ifft(x1,[],dir);
        x1 = fftshift(x1,dir);
    end

    % Decimate k-space
    x1 = reshape(x1, numel(msk), []);
    x1(msk) = 0;
    x1 = reshape(x1, dim);

    % FFT
    for dir=dirs
        x1 = fftshift(x1,dir);
        x1 = fft(x1,[],dir);
    end
    % x1 = fftshift(fftshift(fftshift(x1,1),2),3);

    % Save image
    if size(x,5) == 2
        x1 = cat(5,real(x1),imag(x1));
    end
    x(:,:,:,n,:) = x1;
    
end