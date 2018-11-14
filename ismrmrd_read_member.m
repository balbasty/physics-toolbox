function [obj,nblines] = ismrmrd_read_member(fname, members, idx)
% FORMAT header = ismrmrd_read_member(fname, members, idx)
%
% fname   - Path to the H5 file.
% members - List of members to read.
%           * If a string or cell of strings: names of members to read.
%             One can use dots to nest the selection (e.g. 'head.idx.set')
%           * If a structure: nested description of members to read.
%             (e.g. {'head': {'idx': 'set'}}})
%           * If empty: read all members [default]
% idx     - Indices of elements to read (0-indexing)
%           * If empty: read all elements [default]
%
% Read specific members of a subset k-space lines.

if nargin < 3
    idx = [];
end
if nargin < 2 || isempty(members)
    members = struct;
end

% -------------------------------------------------------------------------
% Build structure describing sub compound type to extract
if ~isstruct(members)
    memberstr = members;
    members   = struct;
    if ~iscell(memberstr)
        memberstr = {memberstr};
    end
    for i=1:numel(memberstr)
        eval(['members.' memberstr{i} ' = struct;']);
    end
end

% -------------------------------------------------------------------------
% Open dataset
file      = H5F.open(fname, 'H5F_ACC_RDONLY','H5P_DEFAULT');
dataset   = H5D.open(file, '/dataset/data');
space     = H5D.get_space(dataset);
type      = H5D.get_type(dataset);

% -------------------------------------------------------------------------
% Build sub compound type
subtype = buildsubtype(type, members);

% -------------------------------------------------------------------------
% Build subspace
if isempty(idx)
    subspace    = 'H5S_ALL';
    memspace    = 'H5S_ALL';
    [~,nblines] = H5S.get_simple_extent_dims(space);
    nblines     = nblines(1);
else
    ndim = H5S.get_simple_extent_ndims(space);
    subspace = space;
    idx      = padarray(idx(:), [0 ndim-1], 0, 'post');
    nblines  = size(idx,1);
    memspace = H5S.create_simple(2, [nblines 1], []);
    H5S.select_elements(subspace,'H5S_SELECT_SET',idx');
end

% -------------------------------------------------------------------------
% Read object
obj = H5D.read(dataset, subtype, memspace, subspace, 'H5P_DEFAULT');

% ==========================================================================
function subtype = buildsubtype(type, members)

membernames = fieldnames(members);

if isempty(membernames)
    subtype = type;
    return
end

subsubtypes   = cell(1,numel(membernames));
size = 0;
for i=1:numel(membernames)
    membername       = membernames{i};
    memberindex      = H5T.get_member_index(type, membername);
    membertype       = H5T.get_member_type(type, memberindex);
    subsubtypes{i}   = buildsubtype(membertype, members.(membername));
    size             = size + H5T.get_size(subsubtypes{i});
end
subtype = H5T.create('H5T_COMPOUND', size);
offset  = 0;
for i=1:numel(membernames)
    H5T.insert(subtype, membernames{i}, offset, subsubtypes{i});
    offset = offset + H5T.get_size(subsubtypes{i});
end