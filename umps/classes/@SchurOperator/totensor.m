function T = totensor(obj)
% TOTENSOR Convert a Schur operator into an MPO tensor.
%
% T = totensor(C)
% Given an operator in Schur form (cell array of (d x d) matrices or
% scalars), returns the corresponding rank-4 tensor, with the same
% convention as MPOs.

d = max(cellfun(@(c) max(size(c)),obj.data(:)));
Cfull = cellfun(@(c) expandterm(c,d),obj.data,'UniformOutput',false);
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
