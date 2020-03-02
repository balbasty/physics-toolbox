function model = init(in, model, opt)
% FORMAT model = b1p.bs.sens.init(in, model, opt)

lat = [size(in.coils,1) size(in.coils,2) size(in.coils,3)];
Nc  = size(in.coils,4);
Np  = size(in.coils,5);

% -------------------------------------------------------------------------
% Process slice-wise to save memory
sumrr = zeros(1,Nc);
sumrx = zeros(1,Nc);
for z=1:lat(3) 
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = single(in.coils(:,:,z,:,:));
    xz = reshape(xz, [], Nc, Np);
    mask = isfinite(xz);
    xz(~mask) = 0;

    % ---------------------------------------------------------------------
    % Load one slice of the transmit field
    bz = single(model.b1(:,:,z));
    bz = reshape(bz, [], 1);
    if opt.b1.log, bz = exp(2*bz); end
    sign = reshape(in.pulse.sign, 1, 1, Np);
    bz = bsxfun(@times, bz, sign);
    bz = bz * in.pulse.factor;
    bz = exp(1i * bz);

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    rz = single(model.mean(:,:,z));
    rz = reshape(rz, [], 1);
    rz = bsxfun(@times, rz, bz);
    clear bz
    
    % ---------------------------------------------------------------------
    % Accumulate
    rx = bsxfun(@times, conj(rz), xz);
    rx(~mask) = 0;
    sumrx = sumrx + sum(sum(rx, 3), 1, 'double');
    clear rx
    
    rr = real(conj(rz).*rz);
    rr = bsxfun(@times, rr, mask);
    sumrr = sumrr + sum(sum(rr, 3), 1, 'double');
    clear rr
end

% -------------------------------------------------------------------------
% Compute ML constant sensitivity
sumrx = sumrx./sumrr;

% -------------------------------------------------------------------------
% Write 
for n=1:Nc
    model.sens(:,:,:,n) = sumrx(n);
end