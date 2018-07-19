function A_new = applyHC_g(C,H,H_left,H_right,varargin)
A_new = ncon({H_left,C,H_right},{[-1,1,3],[1,2],[-2,2,3]});
end