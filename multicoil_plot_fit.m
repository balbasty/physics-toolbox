function multicoil_plot_fit(n, x, s, rho, vs, figname, movie, it)
% FORMAT multicoil_plot_fit(n, x, s, rho, (vs), (figname))

% -------------------------------------------------------------------------
% Set path
path = fileparts(which('multicoil_plot_fit'));
addpath(fullfile(path, 'colormaps'));
addpath(fullfile(path, 'phasemap'));

% -------------------------------------------------------------------------
% Default parameters
if nargin < 8
    it = NaN;
    if nargin < 7
        movie = '';
        if nargin < 6
            figname = sprintf('Multicoil fit');
            if nargin < 5
                vs = [1 1 1];
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
fit = s1.*rho1;
res = x1 - fit;

obsfit    = [x1 fit];
obsfitres = [x1 fit res];


% -------------------------------------------------------------------------
% Magnitude
subplot(2,2,1)
h = imagesc(abs(obsfit));
colormap(h.Parent, viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('magnitude (obs & fit) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Phase
subplot(2,2,2)
h = imagesc(angle(obsfit));
colormap(h.Parent, phasemap(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('phase (obs & fit) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Real
subplot(2,2,3)
h = imagesc(real(obsfitres));
colormap(h.Parent, viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('real (obs & fit) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Real
subplot(2,2,4)
h = imagesc(imag(obsfitres));
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