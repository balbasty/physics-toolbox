%% Load results
load('brain_results_180822.mat');
vs = [8 8 8];

z = ceil(size(x,3)/2);
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
    imagesc(abs(cx(:,:,z,n)));
    caxis([0 maxmagx(n)]);
    colormap(jet(1024))
    daspect(vs);
    axis off
    print(sprintf('brain/magx_%d.png', n), '-dpng');
end

imagesc(abs(crho(:,:,z)));
colormap(jet(1024))
daspect(vs);
axis off
print(sprintf('brain/magrho.png'), '-dpng');
colorbar
print(sprintf('brain/magrho_colorbar.png'), '-dpng');

for n=1:size(cx,4)
    imagesc(abs(cs(:,:,z,n) .* crho(:,:,z)));
    colormap(jet(1024))
    caxis([0 maxmagx(n)]);
    daspect(vs);
    axis off
    colorbar
    print(sprintf('brain/magsrho_%d.png', n), '-dpng');
end

%% Save phase image

for n=1:size(cx,4)
    imagesc(angle(cx(:,:,z,n)));
    caxis([-pi pi]);
    phasemap(1024)
    daspect(vs);
    axis off
    print(sprintf('brain/phasex_%d.png', n), '-dpng');
end

imagesc(angle(crho(:,:,z)));
caxis([-pi pi]);
daspect(vs);
axis off
phasemap(1024)
print(sprintf('brain/phaserho.png'), '-dpng');
colorbar
print(sprintf('brain/phaserho_colorbar.png'), '-dpng');

for n=1:size(cx,4)
    imagesc(angle(cs(:,:,z,n) .* crho(:,:,z)));
    caxis([-pi pi]);
    phasemap(1024)
    daspect(vs);
    axis off
    print(sprintf('brain/phasesrho_%d.png', n), '-dpng');
end

% %% Compute ranges
% 
% res = real(bsxfun(@times,cs,crho) - cx).^2;
% maxrealres = max(res(:));
% res = imag(bsxfun(@times,cs,crho) - cx).^2;
% maximagres = max(res(:));
% clear res
% maxres = max(maxrealres,maximagres);
% 

% %% Save residual
% 
% for n=1:size(cx,4)
%     imagesc(real(cs(:,:,z,n) .* crho(:,:,z) - cx(:,:,z,n)).^2);
%     caxis([0 maxres]);
%     colormap(jet(1024))
%     daspect(vs);
%     axis off
%     print(sprintf('brain/realres_%d.png', n), '-dpng');
%     imagesc(imag(cs(:,:,z,n) .* crho(:,:,z) - cx(:,:,z,n)).^2);
%     caxis([0 maxres]);
%     colormap(jet(1024))
%     daspect(vs);
%     axis off
%     print(sprintf('brain/imagres_%d.png', n), '-dpng');
% end

%% Save sensitivity

for n=1:size(cx,4)
    imagesc(abs(cs(:,:,z,n)));
    colormap(jet(1024))
    daspect(vs);
    axis off
    colorbar
    print(sprintf('brain/mags_%d.png', n), '-dpng');
end

for n=1:size(cx,4)
    imagesc(angle(cs(:,:,z,n)));
    caxis([-pi pi]);
    phasemap(1024)
    daspect(vs);
    axis off
    print(sprintf('brain/phases_%d.png', n), '-dpng');
end