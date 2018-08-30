%% Load results
load('phantom_results_180821_487.mat');
vs = [1 1 1];

z = 1;
N = size(x,4);

%% Convert to complex
crho = double(rho(:,:,:,:,1) + 1i * rho(:,:,:,:,2));
cx   = double(x(:,:,:,:,1) + 1i * x(:,:,:,:,2));
cs   = exp(double(s(:,:,:,:,1) + 1i * s(:,:,:,:,2)));

%% Compute ranges

maxmagx = max(abs(reshape(cx(:,:,z,:), [], N)), [], 1);
maxmags = max(abs(reshape(cs(:,:,z,:), [], N)), [], 1);
minmags = min(abs(reshape(cs(:,:,z,:), [], N)), [], 1);

%% Save magnitude image

for n=1:size(cx,4)
    imagesc(abs(cx(:,:,:,n)));
    caxis([0 maxmagx(n)]);
    colormap(jet(1024))
    daspect(vs);
    axis off
    print(sprintf('images/magx_%d.png', n), '-dpng');
end

imagesc(abs(crho));
colormap(jet(1024))
daspect(vs);
axis off
print(sprintf('images/magrho.png'), '-dpng');
colorbar
print(sprintf('images/magrho_colorbar.png'), '-dpng');

for n=1:size(cx,4)
    imagesc(abs(cs(:,:,:,n) .* crho));
    colormap(jet(1024))
    caxis([0 maxmagx(n)]);
    daspect(vs);
    axis off
    colorbar
    print(sprintf('images/magsrho_%d.png', n), '-dpng');
end

%% Save phase image

for n=1:size(cx,4)
    imagesc(angle(cx(:,:,:,n)));
    caxis([-pi pi]);
    phasemap(1024)
    daspect(vs);
    axis off
    print(sprintf('images/phasex_%d.png', n), '-dpng');
end

imagesc(angle(crho));
caxis([-pi pi]);
daspect(vs);
axis off
phasemap(1024)
print(sprintf('images/phaserho.png'), '-dpng');
colorbar
print(sprintf('images/phaserho_colorbar.png'), '-dpng');

for n=1:size(cx,4)
    imagesc(angle(cs(:,:,:,n) .* crho));
    caxis([-pi pi]);
    phasemap(1024)
    daspect(vs);
    axis off
    print(sprintf('images/phasesrho_%d.png', n), '-dpng');
end

% %% Compute ranges
% 
% res = real(bsxfun(@times,cs,crho) - cx).^2;
% maxrealres = max(res(:));
% res = imag(bsxfun(@times,cs,crho) - cx).^2;
% maximagres = max(res(:));
% clear res
% maxres = 0.5*max(maxrealres,maximagres);
% 
% 
% %% Save residual
% 
% for n=1:size(cx,4)
%     imagesc(real(cs(:,:,:,n) .* crho - cx(:,:,:,n)).^2);
%     caxis([0 maxres]);
%     colormap(jet(1024))
%     daspect(vs);
%     axis off
%     print(sprintf('images/realres_%d.png', n), '-dpng');
%     imagesc(imag(cs(:,:,:,n) .* crho - cx(:,:,:,n)).^2);
%     caxis([0 maxres]);
%     colormap(jet(1024))
%     daspect(vs);
%     axis off
%     print(sprintf('images/imagres_%d.png', n), '-dpng');
% end

%% Save sensitivity

for n=1:size(cx,4)
    imagesc(abs(cs(:,:,:,n)));
    colormap(jet(1024))
    daspect(vs);
    axis off
    colorbar
    print(sprintf('images/mags_%d.png', n), '-dpng');
end

for n=1:size(cx,4)
    imagesc(angle(cs(:,:,:,n)));
    caxis([-pi pi]);
    phasemap(1024)
    daspect(vs);
    axis off
    print(sprintf('images/phases_%d.png', n), '-dpng');
end