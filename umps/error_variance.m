function g = error_variance(A_left,A_right,C,applyHA,applyH2s)
% Approximates the energy variance <(H - E)^2> using the two-site variance
% as explained in arXiv:1711.01104.  
% Nullspaces
N_left = nullspace(A_left,'l');
N_right = nullspace(A_right,'l');
% First projector
A = ncon({A_left,C},{[-1,1,-3],[1,-2]});
B = ncon({conj(N_left),applyHA(A)},{[1,-1,2],[1,-2,2]});
g = trace(B*B');
% Second projector
A2s = ncon({A,A_right},{[-1,1,-2],[1,-4,-3]});
B = ncon({conj(N_left),applyH2s(A2s),conj(N_right)},{[1,-1,2],[1,2,4,3],[-2,3,4]});
g = g + trace(B*B');
end
