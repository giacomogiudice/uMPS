function err = gauge_error(A,C,A_left,A_right)
% Computes the errors |A*C - A_left|^2, |A*C - C*A|^2, |C*A - A_right|^2
% in Frobenius norm
AC = ncon({A_left,C},{[-1,1,-3],[1,-2]});
CA = ncon({C,A_right},{[-1,1],[1,-2,-3]});
err = [norm(AC(:) - A(:))^2,norm(CA(:) - A(:))^2,norm(AC(:) - CA(:))^2];
end