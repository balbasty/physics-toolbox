syms mu_r mu_i c1_r c1_i c2_r c2_i x1_r x1_i x2_r x2_i real
A = sym('A', [2 2], 'real');
%%
com = inline('[r_a -i_a ; i_a r_a]', 'r_a', 'i_a');
x1 = com(x1_r,x1_i); % First complex coil image
x2 = com(x2_r,x2_i); % Second complex coil image
c1 = com(c1_r,c1_i); % First log-bias field
c2 = com(c2_r,c2_i); % Second log-bias field
b1 = expm(c1);       % First exponentiated bias field
b2 = expm(c2);       % Second exponentiated bias field
mu = com(mu_r,mu_i); % Complex mean
%%
r1 = (x1-mu*b1); % First residual
r2 = (x2-mu*b2); % Second residual

% Negative LL = conj(r).' * A * r     (A is symmetric)
L = (r1'*r1)*A(1,1) + (r2'*r2)*A(2,2) + (r1'*r2)*A(1,2) + (r2'*r1)*A(1,2);
L = 0.5*L(1,1); % LL is real

%% Symbolic gradient and Hessian

gc1 = simplify([diff(L, c1_r) diff(L, c1_i)], 1000);
Hc1 = simplify([diff(diff(L, c1_r), c1_r) diff(diff(L, c1_r), c1_i);
                diff(diff(L, c1_i), c1_r) diff(diff(L, c1_i), c1_i)], 1000);
            
%% Rewrite gradient
% 
% gc1 =
%
% /r
%   A1_1*real((mu*b1)*conj(mu*b1 - x1))
% + A1_2*real((mu*b1)*conj(mu*b2 - x2))
%
% => gc1/r = real( (mu*b1) * A(1,:) * [conj(mu*b(:) - x(:))] )
%
% /i
%   A1_1*imag((mu*b1)*conj(x1))
% - A1_2*imag((mu*b1)*conj(mu*b2 - x2))
%
% => gc1/i = - imag( (mu*b1) * A(1,:) * [conj(mu*b(:) - x(:))]  )

% Check
mub1x1 = mu*b1*x1.';
mub1x2 = mu*b1*x2.';
b1b2   = b1*b2.';
gc1r =   A(1,1)*( det(mu)*det(b1)   - mub1x1(1,1) ) ...
       + A(1,2)*( det(mu)*b1b2(1,1) - mub1x2(1,1) );
gc1i = A(1,1)*mub1x1(2,1) + A(1,2)*mub1x2(2,1) - A(1,2)*det(mu)*b1b2(2,1);


assert(isequal(simplify(gc1r-gc1(1)), sym(0))); % OK
assert(isequal(simplify(gc1i-gc1(2)), sym(0))); % OK

%% Rewrite Hessian
% 
% Hc1
%
% /rr
%   A1_1*( 2*abs(mu)^2*abs(b1)^2 - real((mu*b1)*conj(x1)) ) 
% + A1_2*real((mu*b1)*conj(mu*b2 - x2))
%
% => Hc1/rr = real( (mu*b1) * A(1,:) * [conj(mu*b(:) - x(:))] )
%             + A11 * abs(mu)^2 * abs(b1)^2
%
% /ii
%   A1_1*real((mu*b1)*conj(x1))
% + A1_2*real((mu*b1)*conj(x2))
% - A1_2*real((mu*b1)*conj(mu*b2))
%
% => Hc1/ii = - real( (mu*b1) * A(1,:) * [conj(mu*b(:) - x(:))]  )
%               + A11 * abs(mu)^2 * abs(b1)^2
%
% /ir
%   A1_1*imag((mu*b1)*conj(x1)) 
% + A1_2*imag((mu*b1)*conj(x2))
% - A1_2*imag((mu*b1)*conj(mu*b2))
%
% ==> Hc1/ir = gc1/i

% Check
Hc1rr = A(1,1) * (2 * det(mu) * det(b1) - mub1x1(1,1)) ...
      + A(1,2) * (det(mu) * b1b2(1,1) - mub1x2(1,1));
Hc1ii = A(1,1) * mub1x1(1,1) ...
      - A(1,2) * (det(mu) * b1b2(1,1) - mub1x2(1,1));

assert(isequal(simplify(Hc1rr - Hc1(1,1)), sym(0))); % OK
assert(isequal(simplify(Hc1ii - Hc1(2,2)), sym(0))); % OK
assert(isequal(simplify(gc1i  - Hc1(2,1)), sym(0))); % OK

%% Check positivity
%
% Tr(H) = H(1,1) + H(2,2) = 2 * Ajj * abs(mu)^2 * abs(bj)^2 > 0
%
% Det(H) = H(1,1)H(2,2) - H(1,2)^2
%        = (Ajj * abs(mu)^2 * abs(bj)^2)^2 
%          - real( (mu*bj) * A(j,:) * [conj(mu*b(:) - x(:))] )^2
%          - imag( (mu*b1) * A(1,:) * [conj(mu*b(:) - x(:))]  )^2
%        = (Ajj * abs(mu)^2 * abs(bj)^2)^2
%           - abs( (mu*b1) * A(1,:) * [conj(mu*b(:) - x(:))]  )^2

% Recompute gradient/Hessian when mu is optimal
% mu = conj(b)*A*x / (conj(b)*A*b)
%
% => conj(mu*b)*A*(mu*b) = conj(b)*A*x*conj(x)*A*b / (conj(b)*A*b)
%                        = conj(mu*b)*A*x
