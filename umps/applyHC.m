function C_new = applyHC(C,H,B_left,B_right,A_left,A_right,mode)
% APPLYHC Apply a gauging matrix to an effective Hamiltonian.
%
% C_new = applyHC(C,H,B_left,B_right,A_left,A_right,mode)
% Applies the gauging matrix C to the effective Hamiltonian formed by H and
% the environments B_left and B_right. The mode is given by the VUMPS mode.
%
% See also: VUMPS.

if strcmp(mode,'generic') || strcmp(mode,'schur') || strcmp(mode,'multicell')
    C_new = applyHC_generic(C,B_left,B_right);
elseif strcmp(mode,'twosite')
    C_new = applyHC_twosite(C,H,B_left,B_right,A_left,A_right);
else
    error(['Unrecognized mode ' mode '.']);
end
end

function C_new = applyHC_generic(C,B_left,B_right)
C_new = ncon({B_left,C,B_right},{[-1,1,3],[1,2],[-2,2,3]});
end

function C_new = applyHC_twosite(C,H,B_left,B_right,A_left,A_right)
C_new = ncon({A_left,C,A_right,H,conj(A_left),conj(A_right)},{[5,1,3],[1,2],[2,8,4],[6,7,3,4],[5,-1,6],[-2,8,7]});
C_new = C_new + ncon({B_left,C},{[-1,1],[1,-2]});
C_new = C_new + ncon({C,B_right},{[-1,1],[-2,1]});
end
