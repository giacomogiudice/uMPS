function tol_new = update_tol(err,options)
tol_new = min(options.maxtol,max(options.dynamicfactor*err,options.mintol));
end
