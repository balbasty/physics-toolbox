function multicoil_plot_fit(n, x, s, rho, vs, figname)
% FORMAT multicoil_plot_fit(n, x, s, rho, (vs), (figname))

if nargin < 6
    figname = sprintf('Multicoil fit');
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
    rho1 = double(rho(:,:,z,:,:));
    rho1 = rho1(:,:,:,:,1) + 1i * rho1(:,:,:,:,2);
else
    rho1 = double(rho(:,:,z,:));
end
if size(s,5) == 2
    s1 = double(s(:,:,z,n,:));
    s1 = exp(s1(:,:,:,:,1) + 1i * s1(:,:,:,:,2));
else
    s1 = exp(double(s(:,:,z,n)));
end
if size(x,5) == 2
    x1 = double(x(:,:,z,n,:));
    x1 = x1(:,:,:,:,1) + 1i * x1(:,:,:,:,2);
else
    x1 = double(x(:,:,z,n));
end

% -------------------------------------------------------------------------
% Build images
fit = s1.*rho1;
res = x1 - fit;

obsfit = [x1 fit ];
obsfitres = [x1 fit res];


% -------------------------------------------------------------------------
% Magnitude
subplot(2,2,1)
h = imagesc(abs(obsfit));
colormap(h.Parent, viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
title(sprintf('magnitude (obs & fit) %d', n))

% -------------------------------------------------------------------------
% Phase
subplot(2,2,2)
h = imagesc(angle(obsfit));
colormap(h.Parent, phasemap(128));
daspect(h.Parent, vs);
axis off
colorbar
title(sprintf('phase (obs & fit) %d', n))

% -------------------------------------------------------------------------
% Real
subplot(2,2,3)
h = imagesc(real(obsfitres));
colormap(h.Parent, viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
title(sprintf('real (obs & fit & res) %d', n))

% -------------------------------------------------------------------------
% Real
subplot(2,2,4)
h = imagesc(imag(obsfitres));
colormap(h.Parent, viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
title(sprintf('imaginary (obs & fit & res) %d', n))

drawnow