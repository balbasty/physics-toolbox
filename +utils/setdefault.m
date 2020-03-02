function o = setdefault(o, field, value, exists)
% Set default value in option structure. The value is only set if none 
% existed previously.
%
% FORMAT opt = utils.setdefault(opt, field, value)
% opt   - Structure
% field - Hierarchy if fields (cell, or '.' separated)
% value - Default value
%
% EXAMPLE
% >> opt.my.field1 = 1;
% >> opt = setdefault(opt, 'my.field1', 2);
% >> opt.my.field1
% ans = 1
% >> opt = setdefault(opt, 'my.field2', 3);
% >> opt.my.field2
% ans = 3
if nargin < 4, exists = false; end
if ~iscell(field), field = strsplit(field, '.'); end

if isempty(field)
    if ~exists, o = value; end
else
    exists = isfield(o, field{1});
    if ~exists, o.(field{1}) = []; end
    o.(field{1}) = utils.setdefault(o.(field{1}), field(2:end), value, exists);
end