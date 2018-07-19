function T = cell2tensor(C,d)
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