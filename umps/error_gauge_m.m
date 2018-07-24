function err = error_gauge_m(A,C_left,C_right,A_left,A_right,mode)
if nargin == 5
	mode = 'all';
end
AC = ncon({A_left,C_right},{[-1,1,-3],[1,-2]});
CA = ncon({C_left,A_right},{[-1,1],[1,-2,-3]});
if isequal(mode,'all')
	err = max([norm(AC(:) - A(:))^2,norm(CA(:) - A(:))^2,norm(AC(:) - CA(:))^2]);
elseif isequal(mode,'left')
	err = norm(AC(:) - A(:))^2;
elseif isequal(mode,'center')
	err = norm(AC(:) - CA(:))^2;
elseif isequal(mode,'right')
	err = norm(CA(:) - A(:))^2;
end