function X = set_slice(X, dim, ind, S)
% Write a (multidmensional) slice into a (multidimensional) array.
%
% FORMAT array = utils.set_slice(array, dim, ind, slice)
% array - {array}  input array
% dim   - {vector} dimension(s) along which to select a slice
% ind   - {[cell of} vector[s]} indices to select in each dimension
% slice - {array}  slice
%
% EXAMPLES
% >> % Create 5-dimensional array
% >> X = rand(10,10,10,10,10);
% >>
% >> % Create a slice (i.e., a sub-array)
% >> S = rand(2,10,10,1,10);
% >>
% >> % Write the slice into the array. Equivalent to: X([1 2],:,:,5,:) = S
% >> X = utils.set_slice(X, [1 4], {[1 2], 5}, S);
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
X = subsasgn(X,sub,S);
