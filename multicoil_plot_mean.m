function multicoil_plot_mean(rho, C, ll, vs, figname)
% FORMAT multicoil_plot_mean(rho, C, ll, (vs), (figname))

% -------------------------------------------------------------------------
% Set path
path = fileparts(which('multicoil_plot_mean'));
addpath(fullfile(path, 'colormaps'));
addpath(fullfile(path, 'phasemap'));

% -------------------------------------------------------------------------
% Default parameters
if nargin < 5
    figname = 'Multicoil mean';
    if nargin < 4
        vs = [1 1 1];
    end
end

% -------------------------------------------------------------------------
% Find window
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name',figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f); 

% -------------------------------------------------------------------------
% Select slice
z    = ceil(size(rho,3)/2);
rho1 = double(rho(:,:,z,:));

% -------------------------------------------------------------------------
% Magnitude
subplot(2,3,1)
h = imagesc(abs(rho1));
colormap(h.Parent, 'gray')
daspect(h.Parent, vs);
axis off
title('magnitude')

% -------------------------------------------------------------------------
% Phase
subplot(2,3,2)
h = imagesc(angle(rho1));
colormap(h.Parent, phasemap(128));
daspect(h.Parent, vs);
axis off
title('phase')

% -------------------------------------------------------------------------
% Real
subplot(2,3,4)
h = imagesc(real(rho1));
colormap(h.Parent, viridis(128));
daspect(h.Parent, vs);
axis off
title('real')

% -------------------------------------------------------------------------
% Imag
subplot(2,3,5)
h = imagesc(imag(rho1));
colormap(h.Parent, viridis(128));
daspect(h.Parent, vs);
axis off
title('imag')

% -------------------------------------------------------------------------
% Covariance
subplot(2,3,3)
h = imagesc(C); colorbar
colormap(h.Parent, viridis(128));
daspect(h.Parent, [1 1 1]);
title('covariance')

% -------------------------------------------------------------------------
% Log-likelihood
subplot(2,3,6)
cla reset
if size(ll,1) == 2
    yyaxis right
    plot(ll(2,:));
    yyaxis left
end
plot(ll(1,:));
title('log-likelihood')

drawnow