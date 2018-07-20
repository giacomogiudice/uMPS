function [A_left,A_right,schmidt] = canonical(A,D_max)
% Convention is (bond,bond,physical)
opts.issym = false;
opts.isreal = isreal(A);
D = size(A,1);
if nargin == 1
	D_max = D;
else
	D_max = min(D,D_max);
end
% Compute right dominant eigenvector 
[v_r,eta] = eigs(@(v) applyTv(v,A,1,A,'r'),D^2,1,'lm',opts);
v_r = reshape(v_r,[D,D]);
v_r = v_r/trace(v_r);
[W,E] = eig((v_r + v_r').'/2);
X = W*sqrt(E);
% Compute left dominant eigenvector 
[v_l,eta] = eigs(@(v) applyTv(v,A,1,A,'l'),D^2,1,'lm',opts);
v_l = reshape(v_l,[D,D]);
v_l = v_l/trace(v_l);
[W,E] = eig((v_l + v_l')/2);
Y = sqrt(E)*W';
% Do SVD decomposition
[U,S,V] = svd(Y*X,'econ');
% Trim and rescale
U = U(:,1:D_max);
S = S(1:D_max,1:D_max);
S = S/trace(S*S');
V = V(:,1:D_max);
% Compute inverses
schmidt = diag(S);
n = sum(schmidt > eps*D);
S_inv = zeros(D,1);
S_inv = 1./schmidt(1:n);
YX_inv = V*diag(S_inv)*U';
X_inv = YX_inv*Y;
Y_inv = X*YX_inv;
% Regauge A
coeff = 1/sqrt(real(eta));
U_l = V'*X_inv;
U_r = Y_inv*U;
A_left = coeff*ncon({S*U_l,A,U_r},{[-1,1],[1,2,-3],[2,-2]});
A_right = coeff*ncon({U_l,A,U_r*S},{[-1,1],[1,2,-3],[2,-2]});
end
