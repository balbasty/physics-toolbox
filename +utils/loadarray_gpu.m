function X = loadarray_gpu(X, convert)
    if nargin > 1
        X = gpuArray(convert(X()));
    else
        X = gpuArray(X());
    end
end