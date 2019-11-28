syms Bt X
syms a A R1 TR R2 TE 'positive'

% S = A*sin(a)*(1-exp(-R1*TR))/(1-cos(a)*exp(-R1*TR))*exp(-R2*TE);
% S = A*a*(1-exp(-R1*TR))/(1-(1-0.5*a^2)*exp(-R1*TR))*exp(-R2*TE);
% S = A*a*(1-exp(-R1*TR))/(1-exp(-R1*TR))*exp(-R2*TE);
S = A*a*exp(Bt)*R1*TR/(0.5*a^2*exp(2*Bt)+R1*TR)*exp(-R2*TE);
L = 0.5*(S-X)^2;

g_A  = diff(L, A);
g_Bt = diff(L, Bt);
g_R1 = diff(L, R1);
g_R2 = diff(L, R2);

HH_A_A   = diff(g_A, A);
HH_Bt_Bt = diff(g_Bt, Bt);
HH_R1_R1 = diff(g_R1, R1);
HH_R2_R2 = diff(g_R2, R2);
HH_A_Bt  = diff(g_A, Bt);
HH_A_R1  = diff(g_A, R1);
HH_A_R2  = diff(g_A, R2);
HH_Bt_R1 = diff(g_Bt, R1);
HH_Bt_R2 = diff(g_Bt, R2);
HH_R1_R2 = diff(g_R1, R2);

g = [g_A; g_Bt; g_R1; g_R2]

HH = [HH_A_A  HH_A_Bt  HH_A_R1  HH_A_R2; ...
      HH_A_Bt HH_Bt_Bt HH_Bt_R1 HH_Bt_R2; ...
      HH_A_R1 HH_Bt_R1 HH_R1_R1 HH_R1_R2; ...
      HH_A_R2 HH_Bt_R2 HH_R1_R2 HH_R2_R2];
  
H = simplify(subs(HH, X, S), 1000)

D = diag(H);
pretty(simplify(sum(D),1000))
pretty(simplify(sum(D([1 3 4])),1000))

