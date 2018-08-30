function multicoil_plot_mean(rho, C, ll, vs, figname)
% FORMAT multicoil_plot_mean(rho, C, ll, (vs), (figname))

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
z = ceil(size(rho,3)/2);

if size(rho,5) == 2
    rho1 = single(rho(:,:,z,:,:));
    rho1 = rho1(:,:,:,:,1) + 1i * rho1(:,:,:,:,2);
else
    rho1 = single(rho(:,:,z,:));
end

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
colormap(h.Parent, 'hsv')
daspect(h.Parent, vs);
axis off
title('phase')

% -------------------------------------------------------------------------
% Real
subplot(2,3,4)
h = imagesc(real(rho1));
daspect(h.Parent, vs);
axis off
title('real')

% -------------------------------------------------------------------------
% Imag
subplot(2,3,5)
h = imagesc(imag(rho1));
daspect(h.Parent, vs);
axis off
title('imag')

% -------------------------------------------------------------------------
% Covariance
subplot(2,3,3)
h = imagesc(C); colorbar
daspect(h.Parent, [1 1 1]);
title('covariance')

% -------------------------------------------------------------------------
% Log-likelihood
subplot(2,3,6)
plot(ll);
title('log-likelihood')

drawnow