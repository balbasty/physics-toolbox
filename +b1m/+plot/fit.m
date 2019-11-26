function fit(n, x, s, rho, msk, senslog, vs, figname, movie, it)
% FORMAT b1m.plot.fit(n, x, s, rho, (vs), (figname))
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Default parameters
if nargin < 10
    it = NaN;
    if nargin < 9
        movie = '';
        if nargin < 8
            figname = sprintf('Multicoil fit');
            if nargin < 7
                vs = [1 1 1];
                if nargin < 6
                    senslog = false;
                    if nargin < 5
                        msk = 1;
                    end
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
if senslog
    z1  = double(s(:,:,z,n));
    s1  = exp(z1);
else
    s1  = double(s(:,:,z,n));
    z1  = log(s1);
end
x1   = double(x(:,:,z,n));

% -------------------------------------------------------------------------
% Build images
truefit   = s1.*rho1;
foldedfit = b1m.adjoint_forward(truefit,msk);
res       = x1 - foldedfit;

obsfit    = [truefit foldedfit x1];
obsfitres = [foldedfit x1 res];


% -------------------------------------------------------------------------
% Magnitude
p = subplot(3,2,1);
h = imagesc(abs(obsfit));
caxis(p, gather([0 max(abs(obsfit(:)))+eps]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('magnitude (fit | folded | obs) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Phase
p = subplot(3,2,2);
h = imagesc(angle(obsfit));
caxis(p, [-pi pi]);
colormap(h.Parent, utils.color.phasemap(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('phase (fit | folded | obs) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Sensitivity: magnitude
p = subplot(3,4,5);
h = imagesc(exp(real(z1)));
cmax = gather(max(abs(s1(:))));
if cmax > 0
    caxis(p, [0 cmax]);
end
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('magnitude (sens) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Sensitivity: real
p = subplot(3,4,6);
h = imagesc(real(s1));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('real (sens) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Sensitivity: imag
p = subplot(3,4,7);
h = imagesc(imag(s1));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('imag (sens) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Sensitivity: phase
p = subplot(3,4,8);
h = imagesc(imag(z1));
caxis(p, [-pi pi]);
colormap(h.Parent, hsv(1024));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('phase (sens) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Real
p = subplot(3,2,5);
h = imagesc(real(obsfitres));
caxis(p, gather([min(real(obsfitres(:))) max(real(obsfitres(:)))+eps]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('real (folded | obs | res) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Real
p = subplot(3,2,6);
h = imagesc(imag(obsfitres));
caxis(p, gather([min(imag(obsfitres(:))) max(imag(obsfitres(:)))+eps]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('imaginary (folded | obs | res) %d', n);
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