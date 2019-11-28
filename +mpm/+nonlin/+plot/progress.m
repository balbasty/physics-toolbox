function f = mpm_plot_progress(out,ll,figname)

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
    
    
    % ---------------------------------------------------------------------
    % z-slice to plot
    z = ceil(size(out.R1.dat,3)/2);
    
    % ---------------------------------------------------------------------
    % A
    subplot(2,3,1);
    imagesc(exp(out.A.dat(:,:,z)));
    colormap('gray');
    caxis([0 1E5]);
    colorbar;
    axis off
    title('A (Proton density)');
    
    % ---------------------------------------------------------------------
    % R1
    subplot(2,3,2);
    imagesc(exp(out.R1.dat(:,:,z)));
    colormap('gray');
    caxis([0 2]);
    colorbar;
    axis off
    title('R1');
    
    % ---------------------------------------------------------------------
    % MT
    subplot(2,3,3);
    imagesc(1./(1+exp(-out.MT.dat(:,:,z))));
    colormap('gray');
    caxis([0 0.05]);
    colorbar;
    axis off
    title('MT');
    
    % ---------------------------------------------------------------------
    % R2
    subplot(2,3,4);
    imagesc(exp(out.R2.dat(:,:,z)));
    colormap('gray');
    caxis([0 80]);
    colorbar;
    axis off
    title('R2^*');
    
    % ---------------------------------------------------------------------
    % MTV weights
    if isfield(out, 'W')
        subplot(2,3,5);
        imagesc(out.W.dat(:,:,z));
        colormap('gray');
        colorbar;
        axis off
        title('MTV weights');
    end
    
    % ---------------------------------------------------------------------
    % Negative log-likelihood
    subplot(2,3,6);
    plot(ll);
    title('Log-likelihood');

    drawnow
end