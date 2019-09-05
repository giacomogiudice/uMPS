function T = cell2tensor(C,d)
% CELL2TENSOR Convert a Schur operator into an MPO tensor.
%
% T = cell2tensor(C,d)
% Given an operator in Schur form (cell array of (d x d) matrices or
% scalars), returns the corresponding rank-4 tensor, with the same
% convention as MPOs.
%
% See also: VUMPS.

Cfull = cellfun(@(c) expandterm(c,d),C,'UniformOutput',false);
T = cell2mat(Cfull);
end

function m = expandterm(c,d)
if isempty(c)
    m = zeros(d);
elseif isscalar(c)
    m = c*eye(d);
else
    m = c;
end
m = reshape(m,[1,1,d,d]);
end
