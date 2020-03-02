function o = opt(o)
% Set default options.
%
% FORMAT opt_out = b1p.b1s.opt(opt_in)
% opt_in  - User-defined option structure [empty]
% opt_out - All user-defined options are kept untouched, the other are set 
%           to their default value.

% --- Default values
if nargin < 1, o = struct; end
o = utils.setdefault(o, 'nbiter',        10);     % Number of iterations
o = utils.setdefault(o, 'b1.log',        true);   % Log-encoding of the B1 field
o = utils.setdefault(o, 'b1.prec.abs',   0);      % B1: precision absolute
o = utils.setdefault(o, 'b1.prec.mem',   0);      % B1: precision membrane
o = utils.setdefault(o, 'b1.prec.ben',   1E5);    % B1: precision bending
o = utils.setdefault(o, 'mean.prec.abs', 0);      % Mean: precision absolute
o = utils.setdefault(o, 'mean.prec.mem', 1E1);    % Mean: precision membrane
o = utils.setdefault(o, 'mean.prec.ben', 0);      % Mean: precision bending
o = utils.setdefault(o, 'sens.prec.abs', 0);      % Sens: precision absolute
o = utils.setdefault(o, 'sens.prec.mem', 0);    % Sens: precision membrane
o = utils.setdefault(o, 'sens.prec.ben', 1E4);      % Sens: precision bending
o = utils.setdefault(o, 'verbose',       1);      % Verbosity level