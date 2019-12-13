function fit(n, x, s, rho, msk, senslog, vs, figname, movie, it)
% FORMAT b1m.plot.fit(n, coils, sens, meanim, [...])
% n       - Coil index
% coils   - Coil images
% sens    - Sensitivities
% meanim  - Mean image
% msk     - Undersampling mask                          [none]
% senslog - Log-encoded sensitivities                   [false]
% vs      - Voxel size                                  [1 1 1]
% figname - Name of the figure                          ['Multicoil fit']
% movie   - Save figure as a movie in provided file     ['']
% it      - Iteration number                            [NaN]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Default parameters
if nargin < 10 || ~isfinite(it)
    it = NaN;
end
if nargin < 9
    movie = '';
end
if nargin < 8 || isempty(figname)
    figname = sprintf('Multicoil fit');
end
if nargin < 7 || isempty(vs)
    vs = [1 1 1];
end
if nargin < 6 || isempty(senslog) || ~isfinite(senslog)
    senslog = false;
end
if nargin < 5
    msk = [];
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
hasmask = ~(isempty(msk) || (isscalar(msk) && msk));
truefit = s1.*rho1;
if ~hasmask
    res       = x1 - truefit;
    obsfit    = [truefit x1];
    obsfitres = [truefit x1 res];
else
    foldedfit = b1m.adjoint_forward(truefit,msk);
    res       = x1 - foldedfit;
    obsfit    = [truefit foldedfit x1];
    obsfitres = [foldedfit x1 res];
end

iscplx = ~isreal(x);
ncol   = 1 + iscplx;
nrow   = 3;

% -------------------------------------------------------------------------
% Magnitude
p = subplot(nrow,ncol,sub2ind([ncol nrow],1,1));
h = imagesc(abs(obsfit));
caxis(p, gather([0 max(abs(obsfit(:)))+eps]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
if hasmask, strtitle = sprintf('magnitude (fit | folded | obs) %d', n);
else,       strtitle = sprintf('magnitude (fit | obs) %d', n); end
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Phase
if iscplx
p = subplot(nrow,ncol,sub2ind([ncol nrow],2,1));
h = imagesc(angle(obsfit));
caxis(p, [-pi pi]);
colormap(h.Parent, utils.color.phasemap(128));
daspect(h.Parent, vs);
axis off
colorbar
if hasmask, strtitle = sprintf('phase (fit | folded | obs) %d', n);
else,       strtitle = sprintf('phase (fit | obs) %d', n); end
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)
end

% -------------------------------------------------------------------------
% Sensitivity: magnitude
if iscplx
    p = subplot(nrow,2*ncol,sub2ind([2*ncol nrow],1,2));
else
    p = subplot(nrow,ncol,sub2ind([ncol nrow],1,2));
end
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
if iscplx
p = subplot(nrow,2*ncol,sub2ind([2*ncol nrow],2,2));
h = imagesc(real(s1));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('real (sens) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)
end

% -------------------------------------------------------------------------
% Sensitivity: imag
if iscplx
p = subplot(nrow,2*ncol,sub2ind([2*ncol nrow],3,2));
h = imagesc(imag(s1));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('imag (sens) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)
end

% -------------------------------------------------------------------------
% Sensitivity: phase
if iscplx
p = subplot(nrow,2*ncol,sub2ind([2*ncol nrow],4,2));
h = imagesc(imag(z1));
caxis(p, [-pi pi]);
colormap(h.Parent, hsv(1024));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('phase (sens) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)
end

% -------------------------------------------------------------------------
% Real
p = subplot(nrow,ncol,sub2ind([ncol nrow],1,3));
h = imagesc(real(obsfitres));
caxis(p, gather([min(real(obsfitres(:))) max(real(obsfitres(:)))+eps]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('real (fit | obs | res) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)

% -------------------------------------------------------------------------
% Imag
if iscplx
p = subplot(nrow,ncol,sub2ind([ncol nrow],2,3));
h = imagesc(imag(obsfitres));
caxis(p, gather([min(imag(obsfitres(:))) max(imag(obsfitres(:)))+eps]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
colorbar
strtitle = sprintf('imaginary (fit | obs | res) %d', n);
if isfinite(it), strtitle = [strtitle sprintf(' [%02d]', it)]; end
title(strtitle)
end

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