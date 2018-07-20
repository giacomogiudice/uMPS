function [x,flag,relres,iter,resvec] = gmres_(Afun,b,options)
	[x,flag,relres,iter,resvec] = gmres(Afun,b,numel(b),options.tol,options.maxit,[],[],options.v0);
end