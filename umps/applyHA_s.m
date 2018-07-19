function A_new = applyHA_s(A,H,H_left,H_right,varargin)
[D,~,d] = size(A);
chi = size(H,1);
A_new = zeros(D,D,d);

for a = 1:chi
	for b = 1:a
		if ~isempty(H{a,b})
			if isscalar(H{a,a}) & H{a,b} ~= 0
				A_new = A_new + H{a,b}*ncon({H_left(:,:,a),A,H_right(:,:,b)},{[-1,1],[1,2,-3],[-2,2]});
			else
				A_new = A_new + ncon({H_left(:,:,a),A,H{a,b},H_right(:,:,b)},{[-1,1],[1,2,3],[-3,3],[-2,2]});
			end
		end
	end
end
end