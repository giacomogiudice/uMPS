function [A_left,A_right,S,eta] = canonical(A)
% This function generates the left and right canonical forms for a given
% MPS tensor A, satisfying the usual canonical constraints. Additionally 
% one has A_left*S = S*A_right and S diagonal, and corresponds to the 
% Schmidt coefficients with normalization Tr(S^2) = 1. S^2 is also the
% right and left fixed points of A_left and A_right correspondingly.
% Convention is (bond,bond,physical)
opts.issym = false;
opts.isreal = isreal(A);
D = size(A,1);
% Compute right dominant eigenvector 
[rho_right,eta] = eigs(@(v) applyTv(v,A,1,A,'r'),D^2,1,'lm',opts);
rho_right = reshape(rho_right,[D,D]);
rho_right = rho_right/trace(rho_right);
[W,E] = eig((rho_right + rho_right').'/2);
X = W*sqrt(E);
% Compute left dominant eigenvector 
[rho_left,eta] = eigs(@(v) applyTv(v,A,1,A,'l'),D^2,1,'lm',opts);
rho_left = reshape(rho_left,[D,D]);
rho_left = rho_left/trace(rho_left);
[W,E] = eig((rho_left + rho_left')/2);
Y = sqrt(E)*W';
% Do SVD decomposition
[U,S,V] = svd(Y*X,'econ');
% Compute inverses
coeff = 1/sqrt(trace(S*S));
S = coeff*S;
schmidt = diag(S);
n = sum(schmidt > eps*D);
S_inv = zeros(D,1);
S_inv = 1./schmidt(1:n);
S_inv = diag(S_inv);
% Regauge A
coeff = coeff/sqrt(eta);
A_left = coeff*ncon({U'*Y,A,X*V*S_inv},{[-1,1],[1,2,-3],[2,-2]});
A_right = coeff*ncon({S_inv*U'*Y,A,X*V},{[-1,1],[1,2,-3],[2,-2]});
end
