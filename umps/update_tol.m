function tol = update_tol(err,options)
% UPDATE_TOL Update tolerance for the internal solvers given the error.
%
% tol = update_tol(err,options)
% Return a new tolerance tol using the global error err and the options of
% the solver. The formula used is 
% tol = min(maxtol,max(dynamicfactor*err,mintol)).
%
% See also: VUMPS, VUMPS_SETTINGS.

tol = min(options.maxtol,max(options.dynamicfactor*err,options.mintol));
end
