function in = coregister(in)

    % ---------------------------------------------------------------------
    % Is 2D?
    is2d = in{1}.dim(3) == 1;
    if is2d
        in.trf = repmat(eye(4), size(in.mat,3));
        in.mat = in.mat0;
        return
    end

    % ---------------------------------------------------------------------
    % REFERENCE (vol == 1)
    ref = spm_vol(in{1}.echoes{1}.dat.fname);
    in{1}.trf = eye(4);
    in{1}.mat = in{1}.trf\in{1}.mat0;
    ref.mat   = in{1}.mat;
    
    % ---------------------------------------------------------------------
    % MOVING (vol > 1)
    for v=2:numel(in)
        in{v}.mat = in{v}.trf\in{v}.mat0;
        mov = spm_vol(in{v}.echoes{1}.dat.fname);
        for fwhm = [21 14 7]
            mov.mat = in{v}.mat;
            q   = spm_coreg(ref,mov,struct('fwhm',[fwhm fwhm]));
            in{v}.trf = spm_matrix(q(:)');
            in{v}.mat = in{v}.trf\in{v}.mat;
            in{v}.trf = in{v}.mat0/in{v}.mat;
        end
    end
    
end