function C_new = applyHC(C,H,B_left,B_right,A_left,A_right)
% APPLYHC Apply a gauging matrix to an effective Hamiltonian.
%
% C_new = applyHC(C,H,B_left,B_right,A_left,A_right)
% Applies the gauging matrix C to the effective Hamiltonian formed by H and
% the environments B_left and B_right.

switch class(H)
    case {'double','SchurOperator','cell'}
        C_new = applyHC_generic(C,B_left,B_right);
    case 'TwoSiteOperator'
        C_new = applyHC_twosite(C,H,B_left,B_right,A_left,A_right);
    otherwise
        error('Unrecognized class of input operator.');
end
end
