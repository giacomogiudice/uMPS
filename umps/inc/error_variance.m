function g = error_variance(A_left,A_right,C,H,B_left,B_right)
% ERROR_VARIANCE Approximates the error in the variance.
%
% g = error_variance(A_left,A_right,C,H,B_left,B_right)
% Approximates the energy variance <(H - E)^2> using the two-site variance
% as explained in this paper on the <a
% href="matlab:web('https://arxiv.org/abs/1711.01104')">arXiv</a>.
%
% See also: NULLSPACE, VUMPS.
switch class(H)
    case 'TwoSiteOperator'
        g = error_variance_twosite(A_left,C,A_right,H,B_left,B_right);
	    return
    case 'cell'
        B_mid = applyT(B_right,A_right,H{2},A_right,'r');
    otherwise
        H = {H,H};
    B_mid = B_right;
end
% Nullspaces
N_left = nullspace(A_left,'l');
N_right = nullspace(A_right,'r');

AC = ncon({A_left,C},{[-1,1,-3],[1,-2]});
G_left = applyT(B_left,N_left,H{1},AC,'l');
G_right = applyT(B_right,N_right,H{2},A_right,'r');
% First projector

P1 = ncon({G_left,B_mid},{[-1,1,2],[-2,1,2]});
% Second projector
P2 = ncon({G_left,G_right},{[-1,1,2],[-2,1,2]});
g = P1(:)'*P1(:) + P2(:)'*P2(:);
end

function g = error_variance_twosite(A_left,C,A_right,H,B_left,B_right)
% Nullspaces
N_left = nullspace(A_left,'l');
N_right = nullspace(A_right,'r');
% First projector
AC = ncon({A_left,C},{[-1,1,-3],[1,-2]});
A_prime = applyHA(AC,H,B_left,B_right,A_left,A_right);
P1 = ncon({conj(N_left),A_prime},{[1,-1,2],[1,-2,2]});
% Second projector
A2s = ncon({AC,A_right},{[-1,1,-2],[1,-4,-3]});
A2s_prime = applyH2s(A2s,H,B_left,B_right);
P2 = ncon({conj(N_left),A2s_prime,conj(N_right)},{[1,-1,2],[1,2,4,3],[-2,3,4]});
g = P1(:)'*P1(:) + P2(:)'*P2(:);
end
