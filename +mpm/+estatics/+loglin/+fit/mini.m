function out = mini(in)

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
                x(end+1) = log(max(echo.mean, eps('single')));
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
        out.(contrasts{i}) = y(i);
    end
    out.R2s = y(end);
        
end