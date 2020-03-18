function out = to_mpm(estatics, pre, coreg)
% Generate MPM maps from ESTATICS fits using rational approximations of the
% spoiled gradient-echo signal equation.
%
% FORMAT out = mpm.estatics.to_mpm(estatics, [pre], [coreg])
% estatics - Output structure returned by `mpm.estatics.<mode>.fit`
% pre      - Structure of precomputed B1+\- maps (`mpm.io.precomputed`)
% coreg    - Coregister pre-computed B1 maps with estatics maps [true]
% out      - Output structure with fields: A, R1, [MT]

if nargin < 3
    coreg = true;
end
if nargin < 2 || isempty(pre)
    pre = struct;
end

opt.folder = fileparts(estatics.PDw.dat.fname);
opt.fname  = '.nii';
prefix = {'A' 'R1'};
values = [0 0];
if isfield(estatics, 'MTw')
    prefix = [prefix {'MT'}];
    values = [values 0];
end
out = mpm.io.output(prefix, estatics.dim(:), estatics.mat, values', 'single', opt);

threshold_R1 = 2;
threshold_A  = 100000;
threshold_MT = 5;

FA_PDw = estatics.PDw.extras.FA;
TR_PDw = estatics.PDw.extras.TR;
FA_T1w = estatics.T1w.extras.FA;
TR_T1w = estatics.T1w.extras.TR;
FA_MTw = estatics.MTw.extras.FA;
TR_MTw = estatics.MTw.extras.TR;

% -------------------------------------------------------------------------
% Read PDw data
% -------------------------------------------------------------------------
B1m = 0;
if isfield(pre, 'B1m')
    if     isfield(pre.B1m, 'PDw'), B1m = pre.B1m.PDw;
    elseif isfield(pre.B1m, 'all'), B1m = pre.B1m.all;
    else,                           B1m = pre.B1m; end
    if coreg, T = utils.coreg({estatics.PDw.dat.fname estatics.PDw.mat}, ...
                              {B1m.struct.fname B1m.mat});
    else,     T = eye(4); end
    B1m = utils.pull(single(B1m.map()), estatics.dim, (T\B1m.mat)\estatics.PDw.mat);
    B1m = log(max(B1m, eps('single')));
end
PDw = exp(single(estatics.PDw.dat()) - B1m);
clear B1m


B1p = 1;
if isfield(pre, 'B1p')
    if     isfield(pre.B1p, 'PDw'), B1p = pre.B1p.PDw;
    elseif isfield(pre.B1p, 'all'), B1p = pre.B1p.all;
    else,                           B1p = pre.B1p; end
    if coreg, T = utils.coreg({estatics.PDw.dat.fname estatics.PDw.mat}, ...
                              {B1p.struct.fname B1p.mat});
    else,     T = eye(4); end
    if strcmpi(B1p.unit, 'p.u.'), norm = 100; else, norm = 1; end
    B1p = utils.pull(single(B1p.map()), estatics.dim, (T\B1p.mat)\estatics.PDw.mat);
    B1p = max(B1p/norm, eps('single'));
end
FA_PDw = FA_PDw * B1p;
clear B1p

% -------------------------------------------------------------------------
% Read T1w data
% -------------------------------------------------------------------------

B1m = 0;
if isfield(pre, 'B1m')
    fprintf('Read B1m (T1w)\n');
    if     isfield(pre.B1m, 'PDw'), B1m = pre.B1m.T1w;
    elseif isfield(pre.B1m, 'all'), B1m = pre.B1m.all;
    else,                           B1m = pre.B1m; end
    if coreg, T = utils.coreg({estatics.T1w.dat.fname estatics.T1w.mat}, ...
                              {B1m.struct.fname B1m.mat});
    else,     T = eye(4); end
    B1m = utils.pull(single(B1m.map()), estatics.dim, (T\B1m.mat)\estatics.T1w.mat);
    B1m = log(max(B1m, eps('single')));
end
fprintf('Read T1w\n');
T1w = exp(single(estatics.T1w.dat()) - B1m);
clear B1m

B1p = 1;
if isfield(pre, 'B1p')
    fprintf('Read B1p (T1w)\n');
    if     isfield(pre.B1p, 'PDw'), B1p = pre.B1p.T1w;
    elseif isfield(pre.B1p, 'all'), B1p = pre.B1p.all;
    else,                           B1p = pre.B1p; end
    if coreg, T = utils.coreg({estatics.T1w.dat.fname estatics.T1w.mat}, ...
                              {B1p.struct.fname B1p.mat});
    else,     T = eye(4); end
    if strcmpi(B1p.unit, 'p.u.'), norm = 100; else, norm = 1; end
    B1p = utils.pull(single(B1p.map()), estatics.dim, (T\B1p.mat)\estatics.T1w.mat);
    B1p = max(B1p/norm, eps('single'));
end
FA_T1w = FA_T1w .* B1p;
clear B1p

% -------------------------------------------------------------------------
% R1
% -------------------------------------------------------------------------
fprintf('Compute R1\n');
R1 = 0.5 * ( T1w .* (FA_T1w ./ TR_T1w) - PDw .* (FA_PDw ./ TR_PDw) ) ... 
     ./ max( (PDw ./ FA_PDw) - (T1w ./ FA_T1w), eps('single') );
R1(R1 < eps('single')) = eps('single');
R1(R1 > threshold_R1) = threshold_R1;
out.R1.dat(:,:,:) = R1;

% -------------------------------------------------------------------------
% A
% -------------------------------------------------------------------------
fprintf('Compute A\n');
A = ( T1w .* PDw ) .* ( TR_T1w .* (FA_PDw ./ FA_T1w) - TR_PDw .* (FA_T1w ./ FA_PDw) ) ... 
    ./ ( PDw .* (TR_PDw .* FA_PDw) - T1w .* (TR_T1w .* FA_T1w) );
A(A < eps('single')) = eps('single');
A(A > threshold_A) = threshold_A;
out.A.dat(:,:,:) = A;


if isfield(estatics, 'MTw')
    % ---------------------------------------------------------------------
    % Read MTw data
    % ---------------------------------------------------------------------
    clear T1w PDw
    B1m = 0;
    if isfield(pre, 'B1m')
        fprintf('Read B1m (MTw)\n');
        if     isfield(pre.B1m, 'PDw'), B1m = pre.B1m.MTw;
        elseif isfield(pre.B1m, 'all'), B1m = pre.B1m.all;
        else,                           B1m = pre.B1m; end
        if coreg, T = utils.coreg({estatics.MTw.dat.fname estatics.MTw.mat}, ...
                                  {B1m.struct.fname B1m.mat});
        else,     T = eye(4); end
        B1m = utils.pull(single(B1m.map()), estatics.dim, (T\B1m.mat)\estatics.MTw.mat);
        B1m = log(max(B1m, eps('single')));
    end
    fprintf('Read MTw\n');
    MTw = exp(single(estatics.MTw.dat()) - B1m);
    clear B1m

    B1p = 1;
    if isfield(pre, 'B1p')
        fprintf('Read B1p (MTw)\n');
        if     isfield(pre.B1p, 'PDw'), B1p = pre.B1p.MTw;
        elseif isfield(pre.B1p, 'all'), B1p = pre.B1p.all;
        else,                           B1p = pre.B1p; end
        if coreg, T = utils.coreg({estatics.MTw.dat.fname estatics.MTw.mat}, ...
                                  {B1p.struct.fname B1p.mat});
        else,     T = eye(4); end
        if strcmpi(B1p.unit, 'p.u.'), norm = 100; else, norm = 1; end
        B1p = utils.pull(single(B1p.map()), estatics.dim, (T\B1p.mat)\estatics.MTw.mat);
        B1p = max(B1p/norm, eps('single'));
    end
    FA_MTw = FA_MTw * B1p;
    clear B1p

    % ---------------------------------------------------------------------
    % MT
    % ---------------------------------------------------------------------
    fprintf('Compute MT\n');
    MT = (FA_MTw .* A ./ MTw - 1) .* R1 .* TR_MTw - 0.5 .* FA_MTw.^2;
    MT = MT * 100;
    MT(MT < eps('single')) = eps('single');
    MT(MT > threshold_MT) = threshold_MT;
    out.MT.dat(:,:,:) = MT;
end

