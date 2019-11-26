function [sens,meanim] = centre(sens,meanim,reg)
% FORMAT [sens,meanim] = b1m.centre(sens,meanim,[reg])
%
% Centre sensitivity fields by dividing them by the geometric mean of their
% magnitudes.

% We know that, at the optimum, sum{alpha_n*s_n} = 0
% To converge faster, and to avoid bias towards the initial mean
% estimate, we enforce this condition at each iteration.
% This is implemented using the log-sum-exp trick.

% Zero-centering the sensitivities changes the conditional term
% of the log-likelihood. However, in practice, the impact is small, 
% and updating the log-likelihood is costly. Therefore, we keep it as 
% is. This means that the log-likelihood changes slightly and might 
% drop during the first iterations.

Nc = size(sens,4);

if nargin < 3
    reg = 1;
end
reg = utils.pad(reg, [Nc-size(reg,1) 2-size(reg,2)], 'replicate', 'post');

meansensr = 0;
meansensi = 0;
for n=1:Nc
    sensr = real(sens(:,:,:,n));
    sensi = imag(sens(:,:,:,n));
    meansensr = meansensr + reg(n,1)*sensr;
    meansensi = meansensi + reg(n,2)*sensi;
end
meansensr = meansensr/sum(reg(:,1));
meansensi = meansensi/sum(reg(:,2));
meansens  = complex(meansensr, meansensi);
for n=1:Nc
    sens(:,:,:,n) = sens(:,:,:,n) - meansens;
end

meanim = meanim .* exp(meansens);

