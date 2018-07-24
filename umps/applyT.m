function W = applyT(M,A_1,O,A_2,direction)
if iscell(A_1)
	W = applyT_m(M,A_1,O,A_2,direction);
	return
elseif iscell(O)
	W = applyT_s(M,A_1,O,A_2,direction);
	return
end
if direction == 'l'
	if isempty(O)
		W = zeros([size(A_1,1),size(A_2,1)]);
	elseif isscalar(O)
		W = O*ncon({M,conj(A_1),A_2},{[1,2],[1,-1,3],[2,-2,3]});
	elseif ndims(M) == 2
		W = ncon({M,conj(A_1),O,A_2},{[1,4],[1,-1,2],[2,3],[4,-2,3]});
	else
		W = ncon({M,conj(A_1),O,A_2},{[1,3,5],[1,-1,2],[3,-3,2,4],[5,-2,4]});
	end
elseif direction == 'r'
	if isempty(O)
		W = zeros([size(A_1,2),size(A_2,2)]);
	elseif isscalar(O)
		W = O*ncon({M,conj(A_1),A_2},{[1,2],[-1,1,3],[-2,2,3]});
	elseif ndims(M) == 2
		W = ncon({M,conj(A_1),O,A_2},{[1,4],[-1,1,2],[2,3],[-2,4,3]});
	else
		W = ncon({M,conj(A_1),O,A_2},{[1,3,5],[-1,1,2],[-3,3,2,4],[-2,5,4]});
	end
else
	error(['Unrecognized direction' direction '.']);
end
end

function M = applyT_m(M,A_1,O,A_2,direction)
N = length(A_1);
if ~iscell(O)
	O = repmat({O},1,N);
end
if direction == 'l'
	for n = 1:N
		M = applyT(M,A_1{n},O{n},A_2{n},'l');
	end
elseif direction == 'r'
	for n = N:(-1):1
		M = applyT(M,A_1{n},O{n},A_2{n},'r');
	end	
else
	error(['Unrecognized direction' direction]);
end
end

function W = applyT_s(M,A_1,O,A_2,direction)
chi = size(O,1);
W = zeros(size(M));
if direction == 'l'
	for a = chi:-1:1
		for b = chi:-1:a
			if ~isempty(O{b,a})
					W(:,:,a) = W(:,:,a) + applyT(M(:,:,b),A_1,O{b,a},A_2,'l');
			end
		end
	end
elseif direction == 'r'
	for a = 1:chi
		for b = 1:a
			if ~isempty(O{a,b})
					W(:,:,a) = W(:,:,a) + applyT(M(:,:,b),A_1,O{a,b},A_2,'r');
			end
		end
	end
else
	error(['Unrecognized direction' direction '.']);
end
end

