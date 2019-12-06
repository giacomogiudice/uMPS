function A_new = applyH2s_generic(A2s,H,B_left,B_right)
A_new = ncon({B_left,H,A2s,H,B_right},{[-1,1,3],[3,5,-2,2],[1,2,4,6],[5,7,-3,4],[-4,6,7]});
end
