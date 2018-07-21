function [A_left_new,A_right_new,C_new] = increasebond(D_new,A_left,A_right,C,applyHC)
[D,~,d] = size(A_left);
assert(D_new <= D,'Bond dimension can only be increased.');
assert(D_new > 2*D,'Increasing the bond dimension by more than twice is not possible.');
A_left_new = zeros(D_new,D_new,d);
A_right_new = zeros(D_new,D_new,d);
C_new = zeros(D_new,D_new);
% Compute nullspace tensors
N_left = nullspace(A_left,'l');
N_right = nullspace(A_right,'r');
% Compute SVD and truncate
M = applyHC(C);
[U,~,V] = svd(M,'econ');
D_cut = min(D_new - D,size(U,2));
U = U(:,1:D_cut);
V = V(:,1:D_cut);
% Populate the new tensors
for k = 1:d
	A_left_new(1:D,:,k) = [A_left(:,:,k),N_left(:,:,k)*U];
	A_right_new(:,1:D,k) = [A_left(:,:,k);V'*N_right(:,:,k)];
end
C_new(1:D,1:D) = C;
end