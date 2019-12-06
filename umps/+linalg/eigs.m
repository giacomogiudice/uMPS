function varargout = eigs(f,s,varargin)
if isscalar(s)
	s = [s,1];
end
n = prod(s);

options = varargin{end};
if isstruct(options) && isfield(options,'v0') && ~isempty(options.v0)
    varargin{end}.v0 = reshape(options.v0,[n,1]);
end
fapply = @(x) reshape(f(reshape(x,s)),[n,1]);
[varargout{1:nargout}] = eigs(fapply,n,varargin{:});
if nargout > 1
	varargout{1} = reshape(varargout{1},s);
end
end


