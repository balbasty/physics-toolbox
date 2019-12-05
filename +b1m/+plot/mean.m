function mean(rho, A, ll, vs, figname)
% FORMAT b1m.plot.mean(rho, A, ll, (vs), (figname))
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


iscplx = ~isreal(rho);
ncol   = 2 + iscplx;
nrow   = 2;

% -------------------------------------------------------------------------
% Magnitude
p = subplot(nrow,ncol,sub2ind([ncol nrow],1,1));
h = imagesc(abs(rho1)); colorbar
caxis(p, gather([0 max(abs(rho1(:)))+eps]));
colormap(h.Parent, 'gray')
daspect(h.Parent, vs);
axis off
title('magnitude')

% -------------------------------------------------------------------------
% Phase
if iscplx
p = subplot(nrow,ncol,sub2ind([ncol nrow],2,1));
h = imagesc(angle(rho1)); colorbar
caxis(p, [-pi pi]);
colormap(h.Parent, utils.color.phasemap(128));
daspect(h.Parent, vs);
axis off
title('phase')
end

% -------------------------------------------------------------------------
% Real
if iscplx
p = subplot(nrow,ncol,sub2ind([ncol nrow],1,2));
h = imagesc(real(rho1));
caxis(p, gather([min(real(rho1(:))) max(real(rho1(:)))+eps]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
title('real')
end

% -------------------------------------------------------------------------
% Imag
if iscplx
p = subplot(nrow,ncol,sub2ind([ncol nrow],2,2));
h = imagesc(imag(rho1));
caxis(p, gather([min(imag(rho1(:))) max(imag(rho1(:)))+eps]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, vs);
axis off
title('imag')
end

% -------------------------------------------------------------------------
% Covariance
C = utils.invPD(A);
p = subplot(nrow,ncol,sub2ind([ncol nrow],2+iscplx,1));
h = imagesc(C); colorbar
caxis(p, gather([min(C(:)) max(C(:))+eps]));
colormap(h.Parent, utils.color.viridis(128));
daspect(h.Parent, [1 1 1]);
title('covariance')

% -------------------------------------------------------------------------
% Log-likelihood
if ~isempty(ll)
    p = subplot(nrow,ncol,sub2ind([ncol nrow],2+iscplx,2));
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