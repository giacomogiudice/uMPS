function A_new = applyH2s_twosite(A2s,H,B_left,B_right)
A_new = ncon({A2s,H},{[-1,1,2,-4],[-2,-3,1,2]});
A_new = A_new + ncon({B_left,A2s},{[-1,1],[1,-2,-3,-4]});
A_new = A_new + ncon({A2s,B_right},{[-1,-2,-3,1],[1,-4]});
end
