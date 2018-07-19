function w = applyTv(v,A_1,O,A_2,direction)
D_1 = size(A_1,1);
D_2 = size(A_2,1);
if ndims(O) == 4
	chi = size(O,1);
	M = reshape(v,[D_1,D_2,chi]);
	w = applyT(M,A_1,O,A_2,direction);
	w = reshape(w,[D_1*D_2*chi,1]);
else
	M = reshape(v,[D_1,D_2]);
	w = applyT(M,A_1,O,A_2,direction);
	w = reshape(w,[D_1*D_2,1]);
end
end