function varargout = bicgstab_(Afun,b,options)
	[varargout{1:nargout}] = bicgstab(Afun,b,options.tol,options.maxit,[],[],options.v0);
end