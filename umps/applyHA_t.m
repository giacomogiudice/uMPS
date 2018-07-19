function A_new = applyHA_t(A,h,H_left,H_right,A_left,A_right)
A_new = ncon({A_left,A,h,conj(A_left)},{[4,1,2],[1,-2,3],[5,-3,2,3],[4,-1,5]});
A_new = A_new + ncon({A,A_right,h,conj(A_right)},{[-1,1,2],[1,5,3],[-3,4,2,3],[-2,5,4]});
A_new = A_new + ncon({H_left,A},{[-1,1],[1,-2,-3]});
A_new = A_new + ncon({A,H_right},{[-1,1,-3],[-2,1]});
end