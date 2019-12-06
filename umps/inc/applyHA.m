function A_new = applyHA(A,H,B_left,B_right,A_left,A_right)
% APPLYHA Apply an MPS tensor to an effective Hamiltonian.
%
% A_new = applyHA(A,H,B_left,B_right,A_left,A_right)
% Applies the tensor A to the effective Hamiltonian formed by H and the
% environments B_left and B_right.

switch class(H)
    case 'double'
        A_new = applyHA_generic(A,H,B_left,B_right);
    case 'SchurOperator'
        A_new = applyHA_schur(A,H,B_left,B_right);
    case 'TwoSiteOperator'
        A_new = applyHA_twosite(A,H,B_left,B_right,A_left,A_right);
    case 'cell'
        if iscell(H)
            A_new = applyHA_schur(A,H,B_left,B_right);
        else
            A_new = applyHA_generic(A,H,B_left,B_right);
        end
    otherwise
        error('Unrecognized class of input operator.');
end
end
