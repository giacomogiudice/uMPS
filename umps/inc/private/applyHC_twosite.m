function C_new = applyHC_twosite(C,H,B_left,B_right,A_left,A_right)
C_new = ncon({A_left,C,A_right,H,conj(A_left),conj(A_right)},{[5,1,3],[1,2],[2,8,4],[6,7,3,4],[5,-1,6],[-2,8,7]});
C_new = C_new + ncon({B_left,C},{[-1,1],[1,-2]});
C_new = C_new + ncon({C,B_right},{[-1,1],[-2,1]});
end
