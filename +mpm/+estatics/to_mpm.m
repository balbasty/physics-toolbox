function mpm = estatics_to_mpm(out, in, pre, mpmfit)

    if nargin < 3
        pre = struct;
    end
    if nargin < 4
        mpmfit = false;
    end

    outfolder = fileparts(out.PDw.dat.fname);
    if mpmfit
        mpm = prepare_output({'A' 'R2' 'R1' 'MT'}, out.dim, out.mat, outfolder, 'mpm.nii', [0 0 0 0]);
    else
        mpm = prepare_output({'A' 'R1' 'MT'}, out.dim, out.mat, outfolder, 'mpm.nii', [0 0 0]);
    end
    
    threshold_R1 = 2;
    threshold_A  = 100000;
    threshold_MT = 5;
    
    
    % ---------------------------------------------------------------------
    % R2
    % ---------------------------------------------------------------------
    R2 = single(out.R2.dat());
    R2(R2 < eps('single')) = eps('single');
    if mpmfit
        R2 = log(R2);
    end
    mpm.R2.dat(:,:,:) = R2;
    clear  R2
    
    % ---------------------------------------------------------------------
    % Read PDw data
    % ---------------------------------------------------------------------
    
    
    lB1r = 0;
    if isfield(pre, 'B1r')
        i_B1r = find(pre.B1r.vol == 1);
        if isempty(i_B1r)
            i_B1r = find(pre.B1r.vol == Inf);
        end
        if ~isempty(i_B1r)
            lB1r = pull(single(pre.B1r.dat{i_B1r}()), out.dim, pre.B1r.mat(:,:,i_B1r)\out.mat(:,:,1));
            lB1r = log(max(lB1r, eps('single')));
        end
    end
    PDw    = exp(single(out.PDw.dat()) - lB1r);
    clear lB1r
    
    
    B1t = 1;
    if isfield(pre, 'B1t')
        i_B1t = find(pre.B1t.vol == 1);
        if isempty(i_B1t)
            i_B1t = find(pre.B1t.vol == Inf);
        end
        if ~isempty(i_B1t)
            B1t = pull(single(pre.B1t.dat{i_B1t}()), out.dim, pre.B1t.mat(:,:,i_B1t)\out.mat(:,:,1));
            B1t = max(B1t/100, eps('single'));
        end
    end
    FA_PDw = in.FA(1) * B1t;
    clear B1t
    
    TR_PDw = in.TR(1);
    
    % ---------------------------------------------------------------------
    % Read T1w data
    % ---------------------------------------------------------------------
    
    lB1r = 0;
    if isfield(pre, 'B1r')
        i_B1r = find(pre.B1r.vol == 2);
        if isempty(i_B1r)
            i_B1r = find(pre.B1r.vol == Inf);
        end
        if ~isempty(i_B1r)
            lB1r = pull(single(pre.B1r.dat{i_B1r}()), out.dim, pre.B1r.mat(:,:,i_B1r)\out.mat(:,:,1));
            lB1r = log(max(lB1r, eps('single')));
        end
    end
    T1w    = exp(single(out.T1w.dat()) - lB1r);
    clear lB1r
    
    B1t = 1;
    if isfield(pre, 'B1t')
        i_B1t = find(pre.B1t.vol == 2);
        if isempty(i_B1t)
            i_B1t = find(pre.B1t.vol == Inf);
        end
        if ~isempty(i_B1t)
            B1t = pull(single(pre.B1t.dat{i_B1t}()), out.dim, pre.B1t.mat(:,:,i_B1t)\out.mat(:,:,1));
            B1t = max(B1t/100, eps('single'));
        end
    end
    FA_T1w = in.FA(2) .* B1t;
    clear B1t
    
    TR_T1w = in.TR(2);
    
    % ---------------------------------------------------------------------
    % R1
    % ---------------------------------------------------------------------
    R1 = 0.5 * ( T1w .* (FA_T1w ./ TR_T1w) - PDw .* (FA_PDw ./ TR_PDw) ) ... 
         ./ max( (PDw ./ FA_PDw) - (T1w ./ FA_T1w), eps('single') );
    R1(R1 < eps('single')) = eps('single');
    R1(R1 > threshold_R1) = threshold_R1;
    if mpmfit
        R1 = log(R1);
    end
    mpm.R1.dat(:,:,:) = R1;
    
    % ---------------------------------------------------------------------
    % A
    % ---------------------------------------------------------------------
    A = ( T1w .* PDw ) .* ( TR_T1w .* (FA_PDw ./ FA_T1w) - TR_PDw .* (FA_T1w ./ FA_PDw) ) ... 
        ./ ( PDw .* (TR_PDw .* FA_PDw) - T1w .* (TR_T1w .* FA_T1w) );
    A(A < eps('single')) = eps('single');
    A(A > threshold_A) = threshold_A;
    if mpmfit
        A = log(A);
    end
    mpm.A.dat(:,:,:) = A;
    
    
    if numel(in.TR) >= 3
        % -----------------------------------------------------------------
        % Read MTw data
        % -----------------------------------------------------------------
        clear T1w PDw
        lB1r = 0;
        if isfield(pre, 'B1r')
            i_B1r = find(pre.B1r.vol == 3);
            if isempty(i_B1r)
                i_B1r = find(pre.B1r.vol == Inf);
            end
            if ~isempty(i_B1r)
                lB1r = pull(single(pre.B1r.dat{i_B1r}()), out.dim, pre.B1r.mat(:,:,i_B1r)\out.mat(:,:,1));
                lB1r = log(max(lB1r, eps('single')));
            end
        end
        MTw    = exp(single(out.MTw.dat()) - lB1r);
        clear lB1r
        
        B1t = 1;
        if isfield(pre, 'B1t')
            i_B1t = find(pre.B1t.vol == 3);
            if isempty(i_B1t)
                i_B1t = find(pre.B1t.vol == Inf);
            end
            if ~isempty(i_B1t)
                B1t = pull(single(pre.B1t.dat{i_B1t}()), out.dim, pre.B1t.mat(:,:,i_B1t)\out.mat(:,:,1));
                B1t = max(B1t/100, eps('single'));
            end
        end
        FA_MTw = in.FA(3) * B1t;
        clear B1t
        
        TR_MTw = in.TR(3);

        % -----------------------------------------------------------------
        % MT
        % -----------------------------------------------------------------
        MT = (FA_MTw .* A ./ MTw - 1) .* R1 .* TR_MTw - 0.5 .* FA_MTw.^2;
        MT = MT * 100;
        MT(MT < eps('single')) = eps('single');
        MT(MT > threshold_MT) = threshold_MT;
        if mpmfit
            MT = MT/100;
            MT = log(MT./(1-MT));
        end
        mpm.MT.dat(:,:,:) = MT;
    else
        mpm.MT.dat(:,:,:) = 0;
    end
    
    if mpmfit
        mpm.dat = cat(4, mpm.A.dat, mpm.R2.dat, mpm.R1.dat, mpm.MT.dat);
        mpm.A.idx  = 1;
        mpm.R2.idx = 2;
        mpm.R1.idx = 3;
        mpm.MT.idx = 4;
    end
end

