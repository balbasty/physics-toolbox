function ok = isrealarray(X)
    function okk = isrealtype(T)
        okk = numel(T) > 7 || strcmpi(T(1:7),'complex');
    end
    if isa(X, 'file_array')
        ok = all(cellfun(@isrealtype, {X.dtype}));
    else
        ok = isreal(X);
    end
end