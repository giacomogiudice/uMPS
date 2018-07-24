function w = applyTv(v,A_1,O,A_2,direction)

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