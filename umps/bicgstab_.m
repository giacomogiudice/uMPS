function varargout = bicgstab_(Afun,b,options)
if exist('options','var') && ~isempty(options) && ~isempty(fieldnames(options))
	[varargout{1:nargout}] = bicgstab(Afun,b,options.tol,options.maxit,[],[],options.v0);
else
	[varargout{1:nargout}] = bicgstab(Afun,b);
end
end