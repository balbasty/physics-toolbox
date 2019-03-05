function X = expand(X, msk, N, P, val)
X1 = X;
val = cast(val, 'like', X1);
switch val
    case 0
        X = zeros(N, P, 'like', X1);
    case 1
        X = ones(N, P, 'like', X1);
    case Inf
        X = Inf(N, P, 'like', X1);
    otherwise
        if numel(val) == 1
            X = val * ones(N, P, 'like', X1);
        elseif numel(val) == P
            X = repmat(val, [N 1]);
        end
end
X(msk,:) = X1; clear X1