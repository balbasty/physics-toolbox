function out = to_mpm_mini(estatics)
% Generate MPM maps from ESTATICS fits using rational approximations of the
% spoiled gradient-echo signal equation.
%
% FORMAT out = mpm.estatics.to_mpm_mini(estatics)
% estatics - Output structure returned by `mpm.estatics.loglin.fit.mini`

% -------------------------------------------------------------------------
% Copy R2*
% -------------------------------------------------------------------------
out.R2s.dat = estatics.R2s.dat;

% -------------------------------------------------------------------------
% Read PDw data
% -------------------------------------------------------------------------
PDw    = exp(estatics.PDw.dat);
FA_PDw = estatics.PDw.FA;
TR_PDw = estatics.PDw.TR;

% -------------------------------------------------------------------------
% Read T1w data
% -------------------------------------------------------------------------
T1w    = exp(estatics.T1w.dat);
FA_T1w = estatics.T1w.FA;
TR_T1w = estatics.T1w.TR;

% -------------------------------------------------------------------------
% R1
% -------------------------------------------------------------------------
out.R1.dat = 0.5 * ( T1w .* (FA_T1w ./ TR_T1w) - PDw .* (FA_PDw ./ TR_PDw) ) ... 
           ./ max( (PDw ./ FA_PDw) - (T1w ./ FA_T1w), double(eps('single')) );
       
% -------------------------------------------------------------------------
% A
% -------------------------------------------------------------------------

out.A.dat = ( T1w .* PDw ) .* ( TR_T1w .* (FA_PDw ./ FA_T1w) - TR_PDw .* (FA_T1w ./ FA_PDw) ) ... 
    ./ ( PDw .* (TR_PDw .* FA_PDw) - T1w .* (TR_T1w .* FA_T1w) );


if isfield(estatics, 'MTw')
    % ---------------------------------------------------------------------
    % Read MTw data
    % ---------------------------------------------------------------------
    MTw    = exp(estatics.MTw.dat);
    FA_MTw = estatics.MTw.FA;
    TR_MTw = estatics.MTw.TR;
    
    % ---------------------------------------------------------------------
    % MT
    % ---------------------------------------------------------------------
    out.MT.dat = (FA_MTw .* out.A.dat ./ MTw - 1) .* out.R1.dat .* TR_MTw - 0.5 .* FA_MTw.^2;
end

