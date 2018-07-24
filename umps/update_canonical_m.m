function [A_left,A_right] = update_canonical_m(A,C_left,C_right)
[D,~,d] = size(A);
% Reshape A into matrices
A_left = reshape(permute(A,[1,3,2]),[d*D,D]);
A_right = reshape(A,[D,d*D]);

[UC_left,~,VC_left] = svd(C_left,'econ');
[UC_right,~,VC_right] = svd(C_right,'econ');
[UL,~,VL] = svd(A_left,'econ');
[UR,~,VR] = svd(A_right,'econ');
A_left = permute(reshape(UL*VL'*VC_right*UC_right',[D,d,D]),[1,3,2]);
A_right = reshape(VC_left*UC_left'*UR*VR',[D,D,d]);

end