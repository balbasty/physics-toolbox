function pre = precomputed(pre,load)
% Create structure of precomputed auxiliary data.
% 
%   This function tries to understand and organise auxiliary data so that 
%   it can easily be used in conjunction with our fitting functions. The
%   structure can be used either when generating maps using rational 
%   approximations (`mpm.estatics.to_mpm`) or when directly fitting  
%   quantitative maps (`mpm.nonlin.fit`).
%   All input files are memory mapped to save RAM, unless the 'load' 
%   option is specified. Files should be nifti. This function tries to 
%   guess important metadata from filenames or file headers, or reads it 
%   from provided info structures.
% 
% FORMAT out = mpm.io.precomputed(in, ['load'|'map'])
% ------
%
% The input structure `in` should have [optional] fields:
% .(B1p|B1m).[<contrast>]
%   .struct (filename to structural image)
%   .map    (filename to mapped transmit/receive field)
%   .unit   ([p.u.]|none)
%
% The output structure `out` will have fields [if input provided]:
% .(B1p|B1m).[<contrast>]
%   .dim        (volume dimensions)
%   .mat0       (original voxel-to-world matrix)
%   .mat        (registered voxel-to-world matrix)
%   .trf        (transformation matrix)
%   .nii(f)     (nifti file[s]) 
%   .unit       (p.u.|none)
%   .struct     (file_array OR array [if 'load']) Structural image
%   .map        (file_array OR array [if 'load']) Mapped field

% TODO: how should the receive/transmit mpa be called?
% B1m/B1p, Bir/B1t, receive/transmit, smap/b1map ...

% -------------------------------------------------------------------------
% Check load/map option
% -------------------------------------------------------------------------
if nargin < 2
    load = false;
end
if ischar(load)
    load = strcmpi(load, 'load');
end

% -------------------------------------------------------------------------
% Create structure
% -------------------------------------------------------------------------

levels0 = fieldnames(pre);
for i=1:numel(levels0)
    fprintf('%d', i);
    level0  = pre.(levels0{i});
    levels1 = fieldnames(level0);
    if any(strcmpi('struct', levels1)) || any(strcmpi('map', levels1))
        pre.(levels0{i}) = fillvol(level0,load);
    else
        for j=1:numel(levels1)
            fprintf('.');
            level1 = level0.(levels1{j});
            pre.(levels0{i}).(levels1{j}) = fillvol(level1,load);
        end
    end
    fprintf(' ');
end
fprintf('\n');

function vol = fillvol(vol,load)
vol.nii    = nifti;
vol.nii(1) = nifti(vol.struct);
vol.nii(2) = nifti(vol.map);
vol.dim    = [size(vol.nii(1).dat) 1];
vol.dim    = vol.dim(1:3);
if strcmpi(vol.nii(1).mat0_intent, 'Scanner')
    vol.mat0 = vol.nii(1).mat0;
elseif strcmpi(vol.nii(1).mat_intent, 'Scanner')
    vol.mat0 = vol.nii(1).mat;
else
    vol.mat0 = vol.nii(1).mat;
end
vol.trf  = eye(4);
vol.mat  = vol.trf\vol.mat0;
if ~isfield(vol, 'unit')
    vol.unit = 'p.u.';
end
vol.struct = vol.nii(1).dat;
if load, vol.struct = vol.struct(); end
vol.map = vol.nii(2).dat;
if load, vol.map = vol.map(); end