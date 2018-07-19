function A_new = applyHC_t(C,h,H_left,H_right,A_left,A_right)
A_new = ncon({A_left,C,A_right,h,conj(A_left),conj(A_right)},{[5,1,3],[1,2],[2,8,4],[6,7,3,4],[5,-1,6],[-2,8,7]});
A_new = A_new + ncon({H_left,C},{[-1,1],[1,-2]});
A_new = A_new + ncon({C,H_right},{[-1,1],[-2,1]});
end