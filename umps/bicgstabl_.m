function varargout = bicgstabl_(Afun,b,options)
	[varargout{1:nargout}] = bicgstabl(Afun,b,options.tol,options.maxit,[],[],options.v0);
end