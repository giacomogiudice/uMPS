function A_new = applyH2s_t(A2s,h,H_left,H_right)
A_new = ncon({A2s,h},{[-1,1,2,-4],[-2,-3,1,2]});
A_new = A_new + ncon({H_left,A2s},{[-1,1],[1,-2,-3,-4]});
A_new = A_new + ncon({A2s,H_right},{[-1,-2,-3,1],[1,-4]});
end