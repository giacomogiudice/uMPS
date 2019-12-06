function [B_left,B_right,energy] ...
    = update_environments_multicell(A_left,A_right,C_left,C_right,H,B_left,B_right,settings)
if ~isempty(C_left)
    settings.advice.C = C_left;
end
if ~isempty(B_left) && size(B_left,1) == size(C_left,1)
    settings.advice.B = B_left;
end
[B_left,energy_left] = fixedblock(H,A_left,'l',settings);
if ~isempty(C_right)
    settings.advice.C = C_right;
end
if ~isempty(B_right) && size(B_right,1) == size(C_right,1)
    settings.advice.B = B_right;
end
[B_right,energy_right] = fixedblock(shift(H,1),shift(A_right,1),'r',settings);
energy = mean([energy_left,energy_right]);
% Fix normalization in case of generic MPO
if isfloat(H{1})
    B_norm = sqrt(abs(ncon({B_left,conj(C_left),C_left,B_right},{[1,4,3],[1,2],[4,5],[2,5,3]})));
    B_left = B_left/B_norm;
    B_right = B_right/B_norm;
    energy = real(energy);
end
end
