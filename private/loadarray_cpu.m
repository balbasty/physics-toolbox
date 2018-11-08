function X = loadarray_cpu(X, convert)
    if nargin < 2
        convert = @(X) X;
    end
    X = convert(X());
end