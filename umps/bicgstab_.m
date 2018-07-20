function [x,flag,relres,iter,resvec] = bicgstab_(Afun,b,options)
	[x,flag,relres,iter,resvec] = bicgstab(Afun,b,options.tol,options.maxit,[],[],options.v0);
end