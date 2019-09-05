function err = error_gauge(A,A_left,A_right,C_left,C_right)
% ERROR_GAUGE Return error in the gauge between the MPS tensors.
%
% err = error_gauge(A,A_left,A_right,C)
% Returns the maximum between |A - A_left*C|,|C*A_right - A|, and 
% |A_left*C - C*A_right|.
% Each error is computed as 1 - |(overlap)|^2. One of the rank-3 tensors 
% can be empty and the corresponding error will not be computed. 
%
% err = error_gauge(A,A_left,A_right,C_left,C_right)
% Returns the maximum between |A - A_left*C_right|,|C_left*A_right - A|, 
% and |A_left*C_right - C_left*A_right|.
%
% See also: VUMPS.

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

err = max(abs(1 - abs(overlaps).^2));
end
