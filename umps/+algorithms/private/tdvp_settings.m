function settings = tdvp_settings(custom_settings)
% Define default settings
settings = struct;
settings.niter = 100;
settings.timestep = -0.1i;
settings.isreal = false;
settings.verbose = true;
settings.integrator = struct;
settings.integrator.handle = @linalg.expmv;
settings.integrator.options = struct;
settings.integrator.options.tol = 1e-12;
settings.eigsolver = struct;
settings.eigsolver.handle = @linalg.eigs;
settings.eigsolver.mode = 'sr';
settings.eigsolver.options = struct;
settings.eigsolver.options.issym = true;
settings.eigsolver.options.isreal = false;
settings.eigsolver.options.tol = eps;
settings.eigsolver.options.dynamictol = false;
settings.linsolver.handle = @linalg.bicgstab;
settings.linsolver.options = struct;
settings.linsolver.options.tol = eps;
settings.linsolver.options.maxit = 100;
settings.linsolver.options.v0 = [];
settings.linsolver.options.dynamictol = false;
% Merge settings with custom ones
if exist('custom_settings','var') && ~isempty(custom_settings)
	settings = mergestruct(settings,custom_settings);
end
% Update eigsolver settings for real matrices
if settings.isreal
	settings.eigsolver.options.isreal = true;
end
% Do some minimal checks
assert(mod(settings.niter,1) == 0,'%s must be an integer.','settings.maxit');

end
