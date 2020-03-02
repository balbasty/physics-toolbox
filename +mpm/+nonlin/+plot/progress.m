function f = progress(out,ll,scl,means,figname)

    % ---------------------------------------------------------------------
    % Get figure object
    if nargin < 5 || isempty(figname)
        figname = 'MPM fit';
    end
    f = findobj('Type', 'Figure', 'Name', figname);
    if isempty(f)
        f = figure('Name', figname, 'NumberTitle', 'off');
    end
    set(0, 'CurrentFigure', f);   
    clf(f);
    
    if nargin < 4
        means = struct;
    end
    if nargin < 3 || isempty(scl)
        scl = ones(1,numel(ll));
    end
    if ~isempty(scl)
        colours = hsv(max(scl));
    end
    
    ncol = 4;
    if isfield(out, 'U')
        nrow = 3;
    else
        nrow = 2;
    end
    i = 1;
    
    % ---------------------------------------------------------------------
    % z-slice to plot
    z = ceil(2*out.dim(3)/3);
    
    % ---------------------------------------------------------------------
    % PD
    if isfield(out, 'logA')
        subplot(nrow,ncol,i);
        imagesc(exp(out.logA.dat(:,:,z)));
        if isfield(means, 'logA')
            caxis([0 2*exp(means.logA.dat)]);
        elseif isfield(means, 'A')
            caxis([0 2*means.A.dat]);
        end
        colormap('gray');
        colorbar;
        axis off
        title('A');
        i = i+1;
    elseif isfield(out, 'A')
        subplot(nrow,ncol,i);
        imagesc(out.A.dat(:,:,z));
        if isfield(means, 'logA')
            caxis([0 2*exp(means.logA.dat)]);
        elseif isfield(means, 'A')
            caxis([0 2*means.A.dat]);
        end
        colormap('gray');
        colorbar;
        axis off
        title('A');
        i = i+1;
    end
    
    % ---------------------------------------------------------------------
    % T1
    if isfield(out, 'logR1')
        subplot(nrow,ncol,i);
        imagesc(exp(out.logR1.dat(:,:,z)));
        if isfield(means, 'logR1')
            caxis([0 2*exp(means.logR1.dat)]);
        elseif isfield(means, 'R1')
            caxis([0 2*means.R1.dat]);
        else
            caxis([0 2]);
        end
        colormap('gray');
        colorbar;
        axis off
        title('R1');
        i = i+1;
    elseif isfield(out, 'R1')
        subplot(nrow,ncol,i);
        imagesc(out.R1.dat(:,:,z));
        if isfield(means, 'logR1')
            caxis([0 2*exp(means.logR1.dat)]);
        elseif isfield(means, 'R1')
            caxis([0 2*means.R1.dat]);
        else
            caxis([0 2]);
        end
        colormap('gray');
        colorbar;
        axis off
        title('R1');
        i = i+1;
    end
    
    % ---------------------------------------------------------------------
    % MT
    if isfield(out, 'logMT')
        subplot(nrow,ncol,i);
        imagesc(1./(1+exp(-out.logMT.dat(:,:,z))));
        if isfield(means, 'logMT')
            caxis([0 2./(1+exp(-means.logMT.dat))]);
        elseif isfield(means, 'MT')
            caxis([0 2*means.MT.dat]);
        else
            caxis([0 0.05]);
        end
        colormap('gray');
        colorbar;
        axis off
        title('MT');
        i = i+1;
    elseif isfield(out, 'MT')
        subplot(nrow,ncol,i);
        imagesc(out.MT.dat(:,:,z));
        if isfield(means, 'logMT')
            caxis([0 2./(1+exp(-means.logMT.dat))]);
        elseif isfield(means, 'MT')
            caxis([0 2*means.MT.dat]);
        else
            caxis([0 0.05]);
        end
        colormap('gray');
        colorbar;
        axis off
        title('MT');
        i = i+1;
    end
    
    % ---------------------------------------------------------------------
    % R2
    if isfield(out, 'logR2s')
        subplot(nrow,ncol,i);
        imagesc(exp(out.logR2s.dat(:,:,z)));
        if isfield(means, 'logR2s')
            caxis([0 2*exp(means.logR2s.dat)]);
        elseif isfield(means, 'R2s')
            caxis([0 2*means.R2s.dat]);
        else
            caxis([0 80]);
        end
        colormap('gray');
        colorbar;
        axis off
        title('R2^*');
        i = i+1;
    elseif isfield(out, 'R2s')
        subplot(nrow,ncol,i);
        imagesc(out.R2s.dat(:,:,z));
        if isfield(means, 'logR2s')
            caxis([0 2*exp(means.logR2s.dat)]);
        elseif isfield(means, 'R2s')
            caxis([0 2*means.R2s.dat]);
        else
            caxis([0 80]);
        end
        colormap('gray');
        colorbar;
        axis off
        title('R2^*');
        i = i+1;
    end
    
    % ---------------------------------------------------------------------
    % Uncertainty
    if isfield(out, 'U')
        for j=1:size(out.U.dat, 4)
            subplot(nrow,ncol,i);
            imagesc(out.U.dat(:,:,z,j));
            colormap('gray');
            colorbar;
            axis off
            title('Uncertainty');
            i = i+1;
        end
    end
        
    % ---------------------------------------------------------------------
    % MTV weights
    if isfield(out, 'W')
        subplot(nrow,ncol,i);
        imagesc(out.W.dat(:,:,z));
        colormap('gray');
        colorbar;
        axis off
        title('MTV weights');
        i = i+1;
    end
    
    % ---------------------------------------------------------------------
    % Negative log-likelihood
    subplot(nrow,ncol,i);
    plot(ll, 'k');
    hold on
    x = 1:numel(ll);
    for s=unique(scl)
        p = plot(x(scl==s),ll(scl==s));
        p.LineStyle = 'none';
        p.Marker = 's';
        p.MarkerEdgeColor = 'none';
        p.MarkerFaceColor = colours(s,:);
    end
    hold off
    title('Log-likelihood');

    drawnow
    
    f = true;
end