function W = reduceMPO(W,ind_array,iter_max)
% Reduces the dimension of an MPO tensor by removing linear dependencies.
% The convention for the index ordering is (left,right,top,bottom)
if ~exist('ind_array','var') || isempty(ind_array)
	ind_array = [1,2,3,4];
end
if ~exist('iter_max','var') || isempty(iter_max)
	iter_max = 8;
end

for iter = 1:iter_max
	improvement = 0;
	for ind = ind_array
		[W,delta] = reduceDirection(W,ind);
		improvement = improvement + delta;
	end
	if improvement == 0
		break
	end
end
end

function [W,delta] = reduceDirection(W,ind)
% The convention is that the reduction is carried on the last index and
% that the carryover matrix is contracted on the next-to-last index
if ind == 1
	p = [3,4,2,1];
elseif ind == 2
	p = [3,4,1,2];
elseif ind == 3
	p = [1,2,4,3];
elseif ind == 4
	p = [1,2,3,4];
else
	error(['Index ' ind 'exceeds 4.']);
end

Wp = permute(W,p);
[Wn,delta] = reduceLastIndex(Wp);
pt(p) = 1:length(p); % Inverse permutation
W = permute(Wn,pt);
end

function [W,delta] = reduceLastIndex(W)
[dleft,dright,dtop,dbot] = size(W);
M = reshape(W,[dleft*dright*dtop,dbot]);
[R,colind] = rref(M);
dnew = length(colind);
Mnew = M(:,colind);
carryover = R(1:dnew,:);
W = reshape(Mnew,[dleft,dright,dtop,dnew]);
W = ncon({carryover,W},{[-3,1],[-1,-2,1,-4]});
delta = dbot - dnew;
end
