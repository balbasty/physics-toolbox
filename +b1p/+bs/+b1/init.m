function model = init(in, model, opt)
% FORMAT model = b1p.bs.b1.init(in, model, opt)

lat = [size(in.coils,1) size(in.coils,2) size(in.coils,3)];
Nc  = size(in.coils,4);
Np  = size(in.coils,5);

% -------------------------------------------------------------------------
% Process slice-wise to save memory
sumrr = 0;
sumrx = 0;
for z=1:lat(3) 
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = single(in.coils(:,:,z,:,:));
    xz = reshape(xz, [], Nc, Np);
    mask = isfinite(xz);
    xz(~mask) = 0;

    % ---------------------------------------------------------------------
    % Load one slice of the receive field
    sz = single(model.sens(:,:,z,:));
    sz = reshape(sz, [], Nc);

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    rz = single(model.mean(:,:,z));
    rz = reshape(rz, [], 1);
    rz = bsxfun(@times, rz, sz);
    clear sz
    
    % ---------------------------------------------------------------------
    % Accumulate
    rxp = bsxfun(@times, conj(rz), xz(:,:,in.pulse.sign==1));
    rxm = bsxfun(@times, rz, conj(xz(:,:,in.pulse.sign==-1)));
    rx  = cat(3, rxp, rxm); clear rxp rxm
    rx(~mask) = 0;
    sumrx = sumrx + sum(rx(:), 'double');
    clear rx
    
    rr = real(conj(rz).*rz);
    rr = bsxfun(@times, rr, mask);
    sumrr = sumrr + sum(rr(:), 'double');
    clear rr
end

% -------------------------------------------------------------------------
% Compute ML constant sensitivity
sumrx = sumrx./sumrr;
sumrx = sqrt(real(log(sumrx)/(1i*in.pulse.factor)));
if opt.b1.log, sumrx = log(sumrx); end

% -------------------------------------------------------------------------
% Write 
for n=1:Nc
    model.b1(:,:,:) = sumrx;
end