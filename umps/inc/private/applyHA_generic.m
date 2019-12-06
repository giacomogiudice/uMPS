function A_new = applyHA_generic(A,H,B_left,B_right)
A_new = ncon({B_left,A,H,B_right},{[-1,1,2],[1,5,3],[2,4,-3,3],[-2,5,4]});
end
