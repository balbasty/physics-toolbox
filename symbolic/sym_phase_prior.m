% Objective

phi = sym('phi', [2 1], 'real'); % 2 voxels
K   = sym('K', [4 2], 'real');   % 2 differential features

L = sum(1-cos(K*phi));

% Differentiate

g1  = simplify(diff(L, 'phi1'));
g2  = simplify(diff(L, 'phi2'));
H11 = simplify(diff(g1, 'phi1'));
H22 = simplify(diff(g2, 'phi2'));
H21 = simplify(diff(g1, 'phi2'));
H12 = simplify(diff(g2, 'phi1'));
g  = [g1 g2]'
H  = [H11 H21; H12 H22]

% Check

gg = K'*sin(K*phi);
simplify(gg-g)

HH = K'*diag(cos(K*phi))*K;
simplify(HH-H)

% 

c = cos(phi)
s = sin(phi)

gc1 = simplify(diff(L, 'cos(phi1)'))