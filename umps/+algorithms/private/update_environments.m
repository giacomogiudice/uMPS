function [B_left,B_right,energy] ...
    = update_environments(A_left,A_right,C,H,B_left,B_right,settings)
if ~isempty(C)
    settings.advice.C = C;
end
if ~isempty(B_left)
    settings.advice.B = B_left;
end
[B_left,energy_left] = fixedblock(H,A_left,'l',settings);
if ~isempty(B_right)
    settings.advice.B = B_right;
end
[B_right,energy_right] = fixedblock(H,A_right,'r',settings);
energy = mean([energy_left,energy_right]);
% Fix normalization in case of generic MPO
if isfloat(H) && ~isequal(class(H),'TwoSiteOperator')
    B_norm = sqrt(abs(ncon({B_left,conj(C),C,B_right},{[1,4,3],[1,2],[4,5],[2,5,3]})));
    B_left = B_left/B_norm;
    B_right = B_right/B_norm;
    energy = real(energy);
end
end
