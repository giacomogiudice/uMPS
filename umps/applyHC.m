function A_new = applyHC(C,H,B_left,B_right,A_left,A_right,mode)
if strcmp(mode,'generic') || strcmp(mode,'schur') || strcmp(mode,'multicell')
	A_new = applyHC_generic(C,B_left,B_right);
elseif strcmp(mode,'twosite')
	A_new = applyHC_twosite(C,H,B_left,B_right,A_left,A_right);
else
	error(['Unrecognized mode ' mode '.']);
end
end

function A_new = applyHC_generic(C,B_left,B_right)
A_new = ncon({B_left,C,B_right},{[-1,1,3],[1,2],[-2,2,3]});
end

function A_new = applyHC_twosite(C,H,B_left,B_right,A_left,A_right)
A_new = ncon({A_left,C,A_right,H,conj(A_left),conj(A_right)},{[5,1,3],[1,2],[2,8,4],[6,7,3,4],[5,-1,6],[-2,8,7]});
A_new = A_new + ncon({B_left,C},{[-1,1],[1,-2]});
A_new = A_new + ncon({C,B_right},{[-1,1],[-2,1]});
end