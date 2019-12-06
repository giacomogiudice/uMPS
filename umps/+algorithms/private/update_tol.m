function settings = update_tol(settings,err)
% UPDATE_TOL Update dynamic tolerance in settings structure.
if settings.eigsolver.options.dynamictol
    settings.eigsolver.options.tol = update_formula(err,settings.eigsolver.options);
end
if settings.linsolver.options.dynamictol
    settings.linsolver.options.tol = update_formula(err,settings.linsolver.options);
end
end

function tol = update_formula(err,options)
tol = min(options.maxtol,max(options.dynamicfactor*err,options.mintol));
end
