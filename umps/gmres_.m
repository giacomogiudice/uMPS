function varargout = gmres_(Afun,b,options)
if exist('options','var') && ~isempty(options) && isempty(fieldnames(options))
	[varargout{1:nargout}] = gmres(Afun,b,numel(b),options.tol,options.maxit,[],[],options.v0);
else
	[varargout{1:nargout}] = gmres(Afun,b,numel(b));
end
end