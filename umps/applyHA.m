function A_new = applyHA(A,H,B_left,B_right,A_left,A_right,mode)
if strcmp(mode,'generic')
	A_new = applyHA_generic(A,H,B_left,B_right);
elseif strcmp(mode,'schur')
	A_new = applyHA_schur(A,H,B_left,B_right);
elseif strcmp(mode,'twosite')
	A_new = applyHA_twosite(A,H,B_left,B_right,A_left,A_right);
elseif strcmp(mode,'multicell')
	if iscell(H)
		A_new = applyHA_schur(A,H,B_left,B_right);
	else
		A_new = applyHA_generic(A,H,B_left,B_right);
	end
else
	error(['Unrecognized mode ' mode '.']);
end
end

function A_new = applyHA_generic(A,H,B_left,B_right)
A_new = ncon({B_left,A,H,B_right},{[-1,1,2],[1,5,3],[2,4,-3,3],[-2,5,4]});
end

function A_new = applyHA_schur(A,H,B_left,B_right)
[D,~,d] = size(A);
chi = size(H,1);
A_new = zeros(D,D,d);
for a = 1:chi
	for b = 1:a
		if ~isempty(H{a,b})
			if isscalar(H{a,b})
				A_new = A_new + H{a,b}*ncon({B_left(:,:,a),A,B_right(:,:,b)},{[-1,1],[1,2,-3],[-2,2]});
			else
				A_new = A_new + ncon({B_left(:,:,a),A,H{a,b},B_right(:,:,b)},{[-1,1],[1,2,3],[-3,3],[-2,2]});
			end
		end
	end
end
end

function A_new = applyHA_twosite(A,h,B_left,B_right,A_left,A_right)
A_new = ncon({A_left,A,h,conj(A_left)},{[4,1,2],[1,-2,3],[5,-3,2,3],[4,-1,5]});
A_new = A_new + ncon({A,A_right,h,conj(A_right)},{[-1,1,2],[1,5,3],[-3,4,2,3],[-2,5,4]});
A_new = A_new + ncon({B_left,A},{[-1,1],[1,-2,-3]});
A_new = A_new + ncon({A,B_right},{[-1,1,-3],[-2,1]});
end
