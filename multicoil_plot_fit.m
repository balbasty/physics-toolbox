function multicoil_plot_fit(n, x, s, rho, msk, vs, figname, movie, it)
% FORMAT multicoil_plot_fit(n, x, s, rho, (vs), (figname))

% -------------------------------------------------------------------------
% Set path
path = fileparts(which('multicoil_plot_fit'));
addpath(fullfile(path, 'colormaps'));
addpath(fullfile(path, 'phasemap'));

% -------------------------------------------------------------------------
% Default parameters
if nargin < 9
    it = NaN;
    if nargin < 8
        movie = '';
        if nargin < 7
            figname = sprintf('Multicoil fit');
            if nargin < 6
                vs = [1 1 1];
                if nargin < 5
                    msk = 1;
                end
            end
        end
    end
end
if isempty(figname)
    figname = sprintf('Multicoil fit');
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
s1   = double(s(:,:,z,n));
s1   = exp(s1);
x1   = double(x(:,:,z,n));

% -------------------------------------------------------------------------
% Build images
truefit   = s1.*rho1;
foldedfit = multicoil_pushpullwrap(truefit,msk);
res       = x1 - foldedfit;

obsfit    = [x1 foldedfit truefit];
obsfitres = [x1 foldedfit res];


% -------------------------------------------------------------------------
% Magnitude
p = subplot(2,2,1);
h = imagesc(abs(obsfit));
caxis(p, [0 max(abs(obsfit(:)))]);
colormap(h.Parent, viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('magnitude (obs & fit) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Phase
p = subplot(2,2,2);
h = imagesc(angle(obsfit));
caxis(p, [-pi pi]);
colormap(h.Parent, phasemap(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('phase (obs & fit) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Real
p = subplot(2,2,3);
h = imagesc(real(obsfitres));
caxis(p, [min(real(obsfitres(:))) max(real(obsfitres(:)))]);
colormap(h.Parent, viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('real (obs & fit) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Real
p = subplot(2,2,4);
h = imagesc(imag(obsfitres));
caxis(p, [min(imag(obsfitres(:))) max(imag(obsfitres(:)))]);
colormap(h.Parent, viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('imaginary (obs & fit) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

drawnow


% -------------------------------------------------------------------------
% Write movie
if ~isempty(movie)
    frame = getframe(f);
    [imind, cm] = rgb2ind(frame2im(frame), 256);
    framerate = 1;
    if ~exist(movie, 'file')
        imwrite(imind, cm, movie, 'gif', 'Loopcount', inf, 'DelayTime', 1/framerate); 
    else
        imwrite(imind, cm, movie, 'gif', 'WriteMode', 'append', 'DelayTime', 1/framerate); 
    end
end