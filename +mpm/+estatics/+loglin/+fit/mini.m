function out = mini(in)
% Mini ESTATICS loglinear fit based on the mean echo intensity across the
% whole volume.
%
%   This function can be used to get a ballpark estimate of the mean 
%   ESTATICS parameters by fitting the ESTATICS model to single values 
%   per acquired echo.
%
% FORMAT out = mini(in)
% in  - Input structure returned by `mpm.io.input`
% out - Output structure with one field per parameter (e.g. 'T1w', 'R2s')
%       and subfields:
%       .dat - Estimated parameter value
%       .FA  - Flip angle (all parameters but R2s)
%       .TR  - Repetion time (all parameters but R2s)

    % ---------------------------------------------------------------------
    % GET MINI VALUES
    % ---------------------------------------------------------------------
    c = {}; % contrast id
    x = []; % log-value
    t = []; % -TE
    for v=1:numel(in)
        vol = in{v};
        if strcmpi(vol.seq, 'SGE')
            for e=1:numel(vol.echoes)
                echo = vol.echoes{e};
                c{end+1} = vol.type;
                x(end+1) = log(max(echo.mean, double(eps('single'))));
                t(end+1) = -echo.TE;
            end
        end
    end
    c = c(:);
    x = x(:);
    t = t(:);
    
    % ---------------------------------------------------------------------
    % GLM PROJECTION MATRIX
    % ---------------------------------------------------------------------
    contrasts = unique(c);
    idx = cellfun(@(x) find(strcmpi(x, contrasts)), c);
    B = [zeros(numel(t),3) t];
    for i=unique(idx(:)')
        B(idx == i,i) = 1;
    end
    B = pinv(B);
    
    % ---------------------------------------------------------------------
    % LOG-LINEAR FIT
    % ---------------------------------------------------------------------
    y = B * x;
    out = struct;
    for i=unique(idx(:)')
        out.(contrasts{i}).dat = y(i);
    end
    out.R2s.dat = y(end);
    
    % ---------------------------------------------------------------------
    % TR /Flip Angle
    % ---------------------------------------------------------------------
    for v=1:numel(in)
        vol = in{v};
        if strcmpi(vol.seq, 'SGE')
            out.(vol.type).FA = vol.FA;
            out.(vol.type).TR = vol.TR;
        end
    end
        
end