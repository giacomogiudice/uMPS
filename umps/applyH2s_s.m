function A_new = applyH2s_s(A2s,H,H_left,H_right)
W = cell2tensor(H,size(A2s,2));
A_new = applyH2s_g(A2s,W,H_left,H_right);
end