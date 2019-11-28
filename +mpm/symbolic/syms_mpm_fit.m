% -------------------------------------------------------------------------
% Check gradients for MPM nonlinear fit
% -------------------------------------------------------------------------

syms TE TR sina cosa 'positive'   % Fixed parameters
syms lR1 lR2 lA ld 'real'         % Optimised parameters (log)
syms X 'real'                     % Observation

% Precompute useful values
A  = exp(lA);               % Proton density
R1 = exp(lR1);              % 1/T1 decay
R2 = exp(lR2);              % 1/T2 decay
d  = 1 / (1 + exp(-ld));    % MT ratio
md = 1 - d;
e1 = exp(-TR*R1);
e2 = exp(-TE*R2);

% Signal fit and negative log-likelihood
% Only one approximation: TR2 \approx 0 in the MT pulse => e(-R1*TR2) = 1
% For sequences without a MT pulse, we just set d = 0
S = A * sina * md * (1 - e1) / (1 - md * cosa * e1) * e2;
R = S - X;
L = 0.5*R^2;

% Compute gradients automatically
g = [diff(L, lA); ...
     diff(L, lR2); ...
     diff(L, lR1); ...
     diff(L, ld)];
H = [diff(L, lA,  lA) diff(L, lA,  lR2) diff(L, lA,  lR1) diff(L, lA,  ld); ...
     diff(L, lR2, lA) diff(L, lR2, lR2) diff(L, lR2, lR1) diff(L, lR2, ld); ...
     diff(L, lR1, lA) diff(L, lR1, lR2) diff(L, lR1, lR1) diff(L, lR1, ld); ...
     diff(L, ld,  lA) diff(L, ld,  lR2) diff(L, ld,  lR1) diff(L, ld,  ld)];
H = subs(H, X, S);

% Derivative of (1-d) w.r.t. ld
dd = - exp(-ld)./(1 + exp(-ld)).^2;
simplify(dd - diff(md, ld))

% Check that our gradients are correct
gg = [S; ...
      -TE * R2 * S; ...
      -TR * R1 * e1 * S * (md * cosa / (1 - md * cosa .* e1) - 1 ./ (1 - e1)); ...
      dd * S / (md * (1 - md * cosa * e1))];
HH = gg * gg.';
gg = gg * R;

simplify(g-gg, 100)
simplify(H-HH, 100)
