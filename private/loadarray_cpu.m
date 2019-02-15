function X = loadarray_cpu(X, convert)
    if nargin > 1
        X = convert(X());
    end
end