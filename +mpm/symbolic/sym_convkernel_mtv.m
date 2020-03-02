syms lam                % mem
syms v0 v1 v2           % x/y/z (1/v^2)
N  = 3;                 % Nb voxels in each dimension
% w  = sym(ones(N^3,1));  % This is the classical stationary case
w = sym('w', [N N N]);  % This is the new non-stationary case

O  = sym(zeros(N));
OO = kron(kron(O,O),O);
I  = sym(eye(N));
G1 = sym(spdiags(repmat([-1 1],N,1),[ 0 1],N,N)); G1(N,1) =  1; % forward difference
G2 = sym(spdiags(repmat([-1 1],N,1),[-1 0],N,N)); G2(1,N) = -1; % backward difference
G  = {G1 G2};

% Membrane energy
LL = sym(zeros(N^3,N^3));
for i=1:2
    Di = kron(I,kron(I,G{i}))*sqrt(v0); % 1st order / x
    LL = LL + lam*(Di.'*diag(w(:))*Di)/2;
end
for j=1:2
    Dj = kron(I,kron(G{j},I))*sqrt(v1); % 1st order / y
    LL = LL + lam*(Dj.'*diag(w(:))*Dj)/2;
end
for k=1:2
    Dk = kron(G{k},kron(I,I))*sqrt(v2); % 1st order / z
    LL = LL + lam*(Dk.'*diag(w(:))*Dk)/2;
end

% Reshape so that first 3: output voxels, second 3: input voxels
LL = reshape(LL, [N N N N N N]);
% Extract central output voxel
c = ceil(N/2);
LL = reshape(LL(c,c,c,:,:,:), [N N N]);
% The output value depends on all input neighbouring values
% We'll have some dependencies on the weight image 
LL = simplify(LL, 100);

% We have to construct 7 voxel-specific convolution weights.
% Each of these weights is obtained by convolving the weight image

wker = sym(zeros(N,N,N,N,N,N));
for i=1:N
for j=1:N
for k=1:N
    if ~isequal(LL(i,j,k), sym(0))
        for ii=1:N
        for jj=1:N
        for kk=1:N
            if ~isequal(LL(ii,jj,kk), sym(0))
                ww = sprintf('w%d_%d_%d', ii, jj, kk);
                wker(i,j,k,ii,jj,kk) = diff(LL(i,j,k),ww);
            end
        end
        end
        end
    end
end
end
end
wker = simplify(wker, 100);
c    = ceil(N/2);
k111 = reshape(wker(c,c,c,:,:,:), [N N N]);
k110 = reshape(wker(c,c,c-1,:,:,:), [N N N]);
k112 = reshape(wker(c,c,c+1,:,:,:), [N N N]);
k101 = reshape(wker(c,c-1,c,:,:,:), [N N N]);
k121 = reshape(wker(c,c+1,c,:,:,:), [N N N]);
k011 = reshape(wker(c-1,c,c,:,:,:), [N N N]);
k211 = reshape(wker(c+1,c,c,:,:,:), [N N N]);


% Convolution kernel for kernel weights:
%
% k000 (3x3x3)
% ------------
% kk000 = lam0 + 
%         lam1*(v0 + v1 + v2) + 
%         lam2*(4*v0^2 + 2*v0*v1 + 2*v0*v2 + 4*v1^2 + 2*v1*v2 + 4*v2^2)
% kk+00 = kk-00 = lam1*(v0/2) + lam2*(v0*v1 + v0*v2 + v0^2)
% kk0+0 = kk0-0 = lam1*(v1/2) + lam2*(v0*v1 + v1*v2 + v1^2)
% kk00+ = kk00- = lam1*(v2/2) + lam2*(v0*v2 + v1*v2 + v2^2)
% kk++0 = kk--0 = lam2*((v0*v1)/2)
% kk+0+ = kk-0- = lam2*((v0*v2)/2)
% kk0++ = kk0-- = lam2*((v1*v2)/2)
%
% 
% k+00 = sym(k-00) (3x3x3)
% ------------------------
% kk000 = kk+00 = lam1*(-v0/2) + lam2*(-(v0*(4*v0 + 2*v1 + 2*v2))/2)
% kk0+0 = kk0-0 = kk++0 = kk+-0 = lam2*(-(v0*v1)/2)
% kk00+ = kk00- = kk+0+ = kk+0- = lam2*(-(v0*v2)/2)
% kk-00 = 0
% kk-+0 = 0
% kk-0+ = 0
%
% k0+0 = sym(k0-0) (3x3x3)
% ------------------------
% kk000 = k0+0 = lam1*(-v1/2) + lam2*(-(v1*(2*v0 + 4*v1 + 2*v2))/2)
% kk+00 = kk-00 = kk++0 = kk-+0 = lam2*(-(v0*v1)/2)
% kk00+ = kk00- = kk0++ = kk0+- = lam2*(-(v1*v2)/2)
% kk0-0 = 0
% kk+-0 = 0
% kk0-+ = 0
%
% k00+ = sym(k00-) (3x3x3)
% ------------------------
% kk000 = k00+ = lam1*(-v2/2) + lam2*(-(v2*(2*v0 + 2*v1 + 4*v2))/2)
% kk00+ = kk00- = kk++0 = kk-+0 = lam2*(-(v0*v2)/2)
% kk00+ = kk00- = kk0++ = kk0+- = lam2*(-(v1*v2)/2)
% kk00- = 0
% kk+0- = 0
% kk0+- = 0
%
% k*00 = sym(k/00)
% ----------------
% kk+00 = lam2*v0^2
%
% k0*0 = sym(k0/0)
% ----------------
% kk0+0 = lam2*v1^2
%
% k00* = sym(k00/)
% ----------------
% kk00+ = lam2*v2^2
%
% k++0 = sym(k--0) = sym(k+-0) = sym(k-+0)
% ----------------------------------------
% kk000 = kk+00 = kk0+0 = kk++0 = lam2*((v0*v1)/2)