function [nlogl, grad] = nlogl_iidstudent_Wzn(Y, W, zn)
% Evaluate negative log-likelihood -log(p(Y|W,z)) and its gradient
% Model: Y = U*inv(W), U ~ i.i.d. student with d.o.f. vector v = exp(zn)
% v is an N x 1 vector with shock-specific d.o.f.
% Depends: kr.m (by Laurent Sorber)
% Marek Jarocinski, 2023-08

v = exp(zn);
T = size(Y,1);
V = repmat(v(:)', T, 1);
U = Y*W;
Util = (U.^2)./V + 1;
logpU = -(V+1)/2.*log(Util);
const = -0.5*log(v) - betaln(0.5, v/2);
nlogl = -size(Y,1)*log(abs(det(W))) - sum(sum(logpU)) - T*sum(const,"all");

if nargout > 1
    % W
    Wtemp = inv(W)';
    A = -(V+1)./V.*U./Util;
    gradw = T*Wtemp(:) + sum(kr(A',Y'),2);
    % v
    temp1 = -0.5*log(Util);
    temp2 = 0.5*(1+V)./(V.^2).*(U.^2)./Util;
    dc = -0.5.*(v.^-1) - 0.5*psi(0.5*v) + 0.5*psi(0.5*(v+1));
    gradv = sum(temp1+temp2)' + T*dc; %*dv, and dv = exp(x)dx = vdx
    grad = [-gradw; -gradv.*v];
end