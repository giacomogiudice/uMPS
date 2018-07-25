function A_new = applyH2s(A2s,H,B_left,B_right,mode)
if strcmp(mode,'generic')
	A_new = applyH2s_generic(A2s,H,B_left,B_right);
elseif strcmp(mode,'schur') | strcmp(mode,'multicell')
	W = cell2tensor(H,size(A2s,2));
	A_new = applyH2s_generic(A2s,W,B_left,B_right);
elseif strcmp(mode,'twosite')
	A_new = applyH2s_twosite(A2s,H,B_left,B_right);
else
	error(['Unrecognized mode ' mode '.']);
end
end

function A_new = applyH2s_generic(A2s,H,B_left,B_right)
A_new = ncon({B_left,H,A2s,H,B_right},{[-1,1,3],[3,5,-2,2],[1,2,4,6],[5,7,-3,4],[-4,6,7]});
end

function A_new = applyH2s_schur(A2s,H,B_left,B_right)
% TODO
end

function A_new = applyH2s_twosite(A2s,H,B_left,B_right)
A_new = ncon({A2s,H},{[-1,1,2,-4],[-2,-3,1,2]});
A_new = A_new + ncon({B_left,A2s},{[-1,1],[1,-2,-3,-4]});
A_new = A_new + ncon({A2s,B_right},{[-1,-2,-3,1],[1,-4]});
end