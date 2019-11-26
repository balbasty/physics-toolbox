function X = ksub(X, odim)
% FORMAT kdata = utils.ksub(kdata, odim)
%
% Extract the central part of k-space of dimensions `odim`

idim = size(X);
idim = utils.pad(idim, [0 numel(odim)-numel(idim)], 1, 'post');
odim = utils.pad(odim, [0 numel(idim)-numel(odim)], NaN, 'post');
odim(~isfinite(odim)) = idim(~isfinite(odim));

if any(odim > idim(1:numel(odim)))
    error('Output dimensions must be smaller than input dimensions')
end

idx = cell(1, numel(idim));
[idx{:}] = deal(':');

for d=1:numel(odim)
    
    id = idim(d);
    od = odim(d);
    
    if id ~= od
        
        centre = floor(id/2)+1;
        pre    = floor(od/2);
        start  = centre - pre;
        stop   = start + od - 1;
        idx{d} = start:stop;
        
    end 
    
end

S = struct('type', '()', 'subs', {idx});
X = subsref(X, S);