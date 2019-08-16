function err = error_gauge(A,A_left,A_right,C_left,C_right)
% Returns the maximum between the errors |A*C - A_left|, |A*C - C*A|,
% |C*A - A_right|. Each "norm" is computed as 1 - |(overlap)|^2.
if ~exist('C_right','var') || isempty(C_right)
	C_right = C_left;
end
overlaps = ones(1,3);
if isempty(A_left) && isempty(A_right)
	error('A_left and A_right cannot be both empty.');
end
if ~isempty(A_left)
	AC = ncon({A_left,C_right},{[-1,1,-3],[1,-2]});
	if ~isempty(A)
		overlaps(1) = AC(:)'*A(:);
	end
end
if ~isempty(A_right)
	CA = ncon({C_left,A_right},{[-1,1],[1,-2,-3]});
	if ~isempty(A)
		overlaps(2) = CA(:)'*A(:);
	end
end
if ~isempty(A_left) && ~isempty(A_right)
	overlaps(3) = AC(:)'*CA(:);
end
% abs(1 - abs(overlaps).^2)
err = max(abs(1 - abs(overlaps).^2));
end

