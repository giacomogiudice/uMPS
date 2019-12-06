function A_new = applyH2s(A2s,H,B_left,B_right)
% APPLYH2S Apply a 2-site MPS tensor to an effective Hamiltonian.
%
% Note: This is legacy code.
%
% A_new = applyH2s(A2s,H,B_left,B_right)
% Applies the 2site effective Hamiltonian to a 2-site tensor A2s.

switch class(H)
    case 'double'
        A_new = applyH2s_generic(A2s,H,B_left,B_right);
    case {'SchurOperator','cell'}
    	W = totensor(H);
    	A_new = applyH2s_generic(A2s,W,B_left,B_right);
    case 'TwoSiteOperator'
        A_new = applyH2s_twosite(A2s,H,B_left,B_right);
    otherwise
        error('Unrecognized class of input operator.');
end
end

function A_new = applyH2s_schur(A2s,H,B_left,B_right)
% TODO
end

