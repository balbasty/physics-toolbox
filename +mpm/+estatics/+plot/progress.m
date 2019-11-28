function f = progress(out,ll,figname)

    % ---------------------------------------------------------------------
    % Get figure object
    if nargin < 3 || isempty(figname)
        figname = 'MPM fit';
    end
    f = findobj('Type', 'Figure', 'Name', figname);
    if isempty(f)
        f = figure('Name', figname, 'NumberTitle', 'off');
    end
    set(0, 'CurrentFigure', f);   
    clf(f);
    
    
    ncol = 4;
    if isfield(out, 'U')
        nrow = 3;
    else
        nrow = 2;
    end
    i = 1;
    
    % ---------------------------------------------------------------------
    % z-slice to plot
    z = ceil(2*size(out.T1w.dat,3)/3);
    
    % ---------------------------------------------------------------------
    % PD
    if isfield(out, 'PDw')
        subplot(nrow,ncol,i);
        imagesc(exp(out.PDw.dat(:,:,z)));
        colormap('gray');
        colorbar;
        axis off
        title('PDw TE=0');
        i = i+1;
    end
    
    % ---------------------------------------------------------------------
    % T1
    if isfield(out, 'T1w')
        subplot(nrow,ncol,i);
        imagesc(exp(out.T1w.dat(:,:,z)));
        colormap('gray');
        colorbar;
        axis off
        title('T1w TE=0');
        i = i+1;
    end
    
    % ---------------------------------------------------------------------
    % MT
    if isfield(out, 'MTw')
        subplot(nrow,ncol,i);
        imagesc(exp(out.MTw.dat(:,:,z)));
        colormap('gray');
        colorbar;
        axis off
        title('MTw TE=0');
        i = i+1;
    end
    
    % ---------------------------------------------------------------------
    % R2
    if isfield(out, 'R2s')
        subplot(nrow,ncol,i);
        imagesc(out.R2s.dat(:,:,z));
        colormap('gray');
        colorbar;
        axis off
        title('R2^* decay');
        i = i+1;
    end
    
    % ---------------------------------------------------------------------
    % Uncertainty
    if isfield(out, 'U')
        subplot(nrow,ncol,i);
        imagesc(out.U.dat(:,:,z,1));
        colormap('gray');
        colorbar;
        axis off
        title('Uncertainty PDw');
        i = i+1;
        
        subplot(nrow,ncol,i);
        imagesc(out.U.dat(:,:,z,2));
        colormap('gray');
        colorbar;
        axis off
        title('Uncertainty T1w');
        i = i+1;
        
        subplot(nrow,ncol,i);
        imagesc(out.U.dat(:,:,z,3));
        colormap('gray');
        colorbar;
        axis off
        title('Uncertainty MTw');
        i = i+1;
        
        subplot(nrow,ncol,i);
        imagesc(out.U.dat(:,:,z,4));
        colormap('gray');
        colorbar;
        axis off
        title('Uncertainty R2^*');
        i = i+1;
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
    plot(ll, '-s', 'MarkerFaceColor','b');
    title('Log-likelihood');

    drawnow
end