function [in,pre] = coregister(in,pre,type)
% Coregister input volumes (and preocmputed maps) together.
%
%   The first echo of each volume is used for co-registration. The first
%   volume of the input structure (`in{1}`) is taken as the reference.
%   If a `pre` structure is provided, structural scans (e.g., 
%   `pre.B1p.struct`) are used for co-registration.
%
% FORMAT [in,[pre]] = mpm.coregister(in,[pre],[type])
% in   - Input structure (obtained from `mpm.io.input`)
% pre  - Precomputed maps structure (obtained from `mpm.io.precomputed`)
% type - What to coregister:
%       (1) channels    [true] 
%       (2) precomputed [true]
%
% To register the precomputed maps, but not the input channels, use 
% type = [false true].
%
% This function alters the fields 'trf' and 'mat' of each volume (but does
% not alter the headers of the files on disk).

if nargin < 3, type = [true true]; end

% -------------------------------------------------------------------------
% REFERENCE DATA (vol == 1)
ref     = spm_vol(in{1}.echoes{1}.dat.fname);
ref.mat = in{1}.mat;

% -------------------------------------------------------------------------
% REGISTER OTHER CHANNELS
if type(1)
    for v=2:numel(in)
        fprintf('Register %s to %s', in{v}.type, in{1}.type);
        in{v}.trf = utils.coreg(ref, {in{v}.echoes{1}.dat.fname, in{v}.mat});
        in{v}.mat = in{v}.trf\in{v}.mat;
        in{v}.trf = in{v}.mat0/in{v}.mat;
    end
end

% -------------------------------------------------------------------------
% PRECOMPUTED
% This bit is ugly. Level2 encodes contrast-specific maps, which should be
% registered to the appropriate volumes. Maps without a level2 (or whose
% level2 is 'all' or 'default') should be registered with the reference
% volume.
if nargin > 1 && ~isempty(pre) && type(2)
    levels1 = fieldnames(pre);
    for i=1:numel(levels1)
        level1 = pre.(levels1{i});
        if isfield(level1, 'struct')
            fprintf('Register %s to %s', levels1{i}, in{1}.type);
            level1.trf = utils.coreg({in{1}.echoes{1}.dat.fname in{1}.mat}, ...
                                     {level1.struct.fname, level1.mat});
            level1.mat = level1.trf\level1.mat;
            level1.trf = level1.mat0/level1.mat;
        else
            levels2 = fieldnames(level1);
            for j=1:numel(levels2)
                level2 = level1.(levels2{j});
                ref = in{1};
                for v=1:numel(in)
                    if strcmpi(in{v}.type, levels2{j})
                        ref = in{v};
                        break
                    end
                end
                if isfield(level2, 'struct')
                    fprintf('Register %s.%s to %s', levels1{i}, levels2{j}, ref.type);
                    level2.trf = utils.coreg({ref.echoes{1}.dat.fname ref.mat}, ...
                                             {level2.struct.fname, level2.mat});
                    level2.mat = level2.trf\level2.mat;
                    level2.trf = level2.mat0/level2.mat;
                end
            end
        end
    end
end