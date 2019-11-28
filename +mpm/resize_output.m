function out = resize_output(out, dim, fields)
% Resize volumes in the output structure.
%
% FORMAT out = resize_output(out, dim, [fields])
% out    - Output structure (created by `mpm.io.output`)
% dim    - New dimensions (max 3 spatial dimensions)
% fields - Output fields to update [default: all]
%
% Note that if generic fields `dim` and `mat` exist, they are also updated.
    
if nargin < 3
    fields = fieldnames(out);
end

% --- resize each volume
for i=1:numel(fields)
    field = fields{i};
    if isa(out.(field), 'nifti') || ...
            (isstruct(out.(field)) && isfield(out.(field), 'dat'))

        dat = out.(field).dat();
        [dat, mat] = utils.upsample(dat, dim);
        out.(field).mat = out.(field).mat * mat;
        if isa(out.(field).dat, 'file_array')
            out.(field).dat.dim(1:3) = dim;
        end
        out.(field).dat(:) = dat(:);
        clear dat
    end
    if isa(out.(field), 'nifti')
        create(out.(field));
    end
end

% --- shared fields
if isfield(out, 'dim')
    dim0    = out.dim;
    out.dim = dim;
    if isfield(out, 'mat')
        scale = dim0 ./ dim;
        mat = [diag(scale) 0.5*(1-scale(:));
               zeros(1,3) 1];
        out.mat = out.mat * mat;
    end
end


