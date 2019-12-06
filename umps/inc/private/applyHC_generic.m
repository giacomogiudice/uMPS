function C_new = applyHC_generic(C,B_left,B_right)
C_new = ncon({B_left,C,B_right},{[-1,1,3],[1,2],[-2,2,3]});
end
