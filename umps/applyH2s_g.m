function A_new = applyH2s_g(A2s,H,H_left,H_right)
A_new = ncon({H_left,H,A2s,H,H_right},{[-1,1,3],[3,5,-2,2],[1,2,4,6],[5,7,-3,4],[-4,6,7]});
end