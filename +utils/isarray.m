function ok = isarray(X)
% Check that a variable is a valid array (numeric OR logical OR file_array)
%
% FORMAT ok = utils.isarray(x)
    ok = isnumeric(X) || islogical(X) || isa(X, 'file_array');
end