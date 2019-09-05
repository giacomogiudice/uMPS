function w = applyTv(v,A_1,O,A_2,direction)
% APPLYTV Wrapper for applyT, but accepts a vector and returns a vector.
%
% w = applyTv(v,A_1,O,A_2,direction) v is reshaped into a matrix
% (column-major) and the transfer matrix is applied. The output vector w is
% the resulting matrix reshaped into a vector.
%
% See also APPLYT.


if iscell(A_1)
    D_1 = size(A_1{1},1);
    D_2 = size(A_2{1},1);
else
    D_1 = size(A_1,1);
    D_2 = size(A_2,1);
end
M = reshape(v,D_1,D_2,[]);
M = applyT(M,A_1,O,A_2,direction);
w = reshape(M,[],1);

end
