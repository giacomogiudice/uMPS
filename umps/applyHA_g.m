function A_new = applyHA_g(A,H,H_left,H_right,varargin)
A_new = ncon({H_left,A,W,H_right},{[-1,1,2],[1,5,3],[2,4,-3,3],[-2,5,4]});
end