function S = select_slice(X, dim, ind)
% Select a (multidimensional) slice from a (multidimensional) array.
%
% FORMAT slice = utils.select_slice(array, dim, ind)
% array - {array}  input array
% dim   - {vector} dimension(s) along which to select a slice
% ind   - {[cell of] vector[s]} indices to select in each dimension
% slice - {array}  slice
%
% EXAMPLES
% >> % Create 5-dimensional array
% >> X = rand(10,10,10,10,10);
% >>
% >> % Select the 2nd slice along the 3rd dimension
% >> % Equivalent to: S = X(:,:,2,:,:)
% >> S = utils.select_slice(X, 3, 2);
% >>
% >> % Select the 1st and 2nd slice along the 1st dimension and the 5th
% >> % along the 4th dimension
% >> % Equivalent to: S = X([1 2],:,:,5,:)
% >> S = utils.select_slice(X, [1 4], {[1 2], 5});
sizeX    = size(X);
sub      = struct;
sub.type = '()';
sub.subs = repmat({':'}, [1 numel(sizeX)]);
if iscell(ind)
    assert(numel(dim) == numel(ind));
    for d=1:numel(dim)
        sub.subs{dim(d)} = ind{d};
    end
else
    sub.subs{dim} = ind;
end
S = subsref(X,sub);
