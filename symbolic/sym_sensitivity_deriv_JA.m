% Define some real values
syms mu_r mu_i b_r b_i f1_r f1_i f2_r f2_i real

% Represent complex numbers by 2x2 matrices
com = inline('[r_a -i_a; i_a r_a]','r_a','i_a');
f1 = com(f1_r,f1_i);
f2 = com(f2_r,f2_i);
b  = com(b_r,b_i);
mu = com(mu_r,mu_i);

% Compute objective function
r1 = (f1-mu*expm( b));
r2 = (f2-mu*expm(-b));
% This should be real, so a scaled identity matrix
E  = (1/2)*(r1'*r1) + (1/2)*(r2'*r2);
E  = E(1,1); % Represent more simply

% 1st derivatives w.r.t. mu
gmu = simplify([diff(E, mu_r)
                diff(E, mu_i)],1000)

% Hessian w.r.t. mu
Hmu = simplify([diff(diff(E,mu_r),mu_r) diff(diff(E,mu_r),mu_i)
                diff(diff(E,mu_i),mu_r) diff(diff(E,mu_i),mu_i)],1000)

% 1st derivatives w.r.t b
gb  = simplify([diff(E, b_r)
                diff(E, b_i)],1000)

% Hessian w.r.t. b - I have no idea if it is positive definite
Hb  = simplify([diff(diff(E,b_r),b_r) diff(diff(E,b_r),b_i)
                diff(diff(E,b_i),b_r) diff(diff(E,b_i),b_i)],1000)