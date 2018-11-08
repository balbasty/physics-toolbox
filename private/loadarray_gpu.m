function X = loadarray_gpu(X, convert)
    if nargin < 2
        convert = @(X) X;
    end
    X = gpuArray(convert(X()));
end