function varargout = gmres_(Afun,b,options)
	[varargout{1:nargout}] = gmres(Afun,b,numel(b),options.tol,options.maxit,[],[],options.v0);
end