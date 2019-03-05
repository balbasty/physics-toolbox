function mean(rho, C, ll, vs, figname)
% FORMAT b1m.plot.mean(rho, C, ll, (vs), (figname))
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

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
p = subplot(2,3,1);
h = imagesc(abs(rho1)); colorbar
caxis(p, gather([0 max(abs(rho1(:)))]));
colormap(h.Parent, 'gray')
daspect(h.Parent, vs);
axis off
title('magnitude')

% -------------------------------------------------------------------------
% Phase
p = subplot(2,3,2);
h = imagesc(angle(rho1)); colorbar
caxis(p, [-pi pi]);
colormap(h.Parent, utils.color.phasemap(128));
daspect(h.Parent, vs);
axis off
title('phase')

% -------------------------------------------------------------------------
% Real
p = subplot(2,3,4);
h = imagesc(real(rho1));
caxis(p, gather([min(real(rho1(:))) max(real(rho1(:)))]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
title('real')

% -------------------------------------------------------------------------
% Imag
p = subplot(2,3,5);
h = imagesc(imag(rho1));
caxis(p, gather([min(imag(rho1(:))) max(imag(rho1(:)))]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
title('imag')

% -------------------------------------------------------------------------
% Covariance
p = subplot(2,3,3);
h = imagesc(C); colorbar
caxis(p, gather([min(C(:)) max(C(:))]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, [1 1 1]);
title('covariance')

% -------------------------------------------------------------------------
% Log-likelihood
if ~isempty(ll)
    subplot(2,3,6)
    cla reset
    if size(ll,1) == 2
        yyaxis right
        plot(ll(2,:));
        yyaxis left
    end
    plot(ll(1,:));
    title('log-likelihood')
end

drawnow