function [A_left,A_right] = update_canonical(A,C)
[D,~,d] = size(A);
% Reshape A into matrices
A_left = reshape(permute(A,[1,3,2]),[d*D,D]);
A_right = reshape(A,[D,d*D]);

method = 'svd';
if isequal(method,'qr')
	% Perform QR decomposition for new A_left
	[QC,~] = qrpos(C);
	[QL,~] = qrpos(A_left);
	A_left = permute(reshape(QL*QC',[D,d,D]),[1,3,2]);
	% Perform LQ decomposition for new A_right
	[~,QC] = lqpos(C);
	[~,QR] = lqpos(A_right);
	A_right = reshape(QC'*QR,[D,D,d]);
elseif isequal(method,'svd')
	[UC,S,VC] = svd(C,'econ');
	[UL,~,VL] = svd(A_left,'econ');
	[UR,~,VR] = svd(A_right,'econ');
	A_left = permute(reshape(UL*VL'*VC*UC',[D,d,D]),[1,3,2]);
	A_right = reshape(VC*UC'*UR*VR',[D,D,d]);
end
end

function [Q,R] = qrpos(A)
% Computes the QR decomposition with R having a positive diagonal
[Q,R] = qr(A,0);
U = diag(sign(diag(R)));
Q = Q*U;
R = U*R;
end

function [L,Q] = lqpos(A)
% Computes the LQ decomposition with L having a positive diagonal
[Q,L] = qr(A.',0);
U = diag(sign(diag(L)));
L = L.'*U;
Q = U*Q.';
end

function [U,P,P_prime] = polard(M)
% Computes the left and right polar decomposition of M as
% M = U*P = P_prime*U
[W,S,V] = svd(M,'econ');

U = W*V';
P = V*S*V';
P = (P + P')/2;
P_prime = W*S*W';
P_prime = (P_prime + P_prime')/2;
end