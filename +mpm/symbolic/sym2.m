syms X A
syms R2 TE 'positive'

S = exp(A-R2*TE);
L = 0.5*(S-X)^2;

g_A  = diff(L, A);
g_R2 = diff(L, R2);

HH_A_A   = diff(g_A, A);
HH_R2_R2 = diff(g_R2, R2);
HH_A_R2  = diff(g_A, R2);

g = [g_A; g_R2]

HH = [HH_A_A  HH_A_R2; ...
      HH_A_R2  HH_R2_R2];
  
H = simplify(subs(HH, X, S), 1000)

%%

syms X A
syms R2 TE 'positive'

S = A-R2*TE;
L = 0.5*(S-X)^2;

g_A  = diff(L, A);
g_R2 = diff(L, R2);

HH_A_A   = diff(g_A, A);
HH_R2_R2 = diff(g_R2, R2);
HH_A_R2  = diff(g_A, R2);

g = [g_A; g_R2]

HH = [HH_A_A  HH_A_R2; ...
      HH_A_R2  HH_R2_R2];
  
H = simplify(subs(HH, X, S), 1000)
