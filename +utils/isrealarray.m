function ok = isrealarray(X)
% Check that a variable is an array of real number
% (numeric OR logical OR file_array)
%
% FORMAT ok = utils.isrealarray(x)
    function okk = isrealtype(T)
        okk = numel(T) > 7 || strcmpi(T(1:7),'complex');
    end
    if isa(X, 'file_array')
        ok = all(cellfun(@isrealtype, {X.dtype}));
    else
        ok = utils.isarray(X) && isreal(X);
    end
end