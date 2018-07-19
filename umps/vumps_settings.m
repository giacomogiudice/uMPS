function settings = vumps_settings(custom_settings)
% Define default settings
settings = struct;
settings.mode = 'schur';
settings.maxit = 100;
settings.tol = 1e-12;
settings.isreal = false;
settings.verbose = true;
settings.initial = struct;
settings.eigsolver = struct;
settings.eigsolver.handle = @eigs;
settings.eigsolver.mode = 'sr';
settings.eigsolver.options = struct;
settings.eigsolver.options.issym = true;
settings.eigsolver.options.isreal = false;
settings.eigsolver.options.tol = eps;
settings.eigsolver.options.dynamictol = false;
settings.eigsolver.options.dynamicfactor = 1e-2;
settings.eigsolver.options.mintol = eps;
settings.eigsolver.options.maxtol = 1e-3;
settings.linsolver.handle = @gmres;
settings.linsolver.options.tol = eps;
settings.linsolver.options.maxit = 500;
settings.linsolver.options.dynamictol = false;
settings.linsolver.options.dynamicfactor = 1e-2;
settings.linsolver.options.mintol = eps;
settings.linsolver.options.maxtol = 1e-3;
% Merge settings with custom ones
if nargin == 1
	settings = mergestruct(settings,custom_settings);
end
% Update eigsolver settings for real matrices
if settings.isreal 
	settings.eigsolver.options.isreal = true;
	settings.eigsolver.mode = 'sa';
end
% Do some minimal checks
assert(mod(settings.maxit,1) == 0);
assert(settings.tol <= 1 & settings.tol >= 0);
assert(settings.eigsolver.options.dynamicfactor <= 1 & settings.eigsolver.options.dynamicfactor >= 0);
assert(settings.eigsolver.options.maxtol >= settings.eigsolver.options.mintol);
assert(mod(settings.linsolver.options.maxit,1) == 0);
assert(settings.linsolver.options.tol <= 1 & settings.linsolver.options.tol >= 0);
assert(settings.linsolver.options.dynamicfactor <= 1 & settings.linsolver.options.dynamicfactor >= 0);
assert(settings.linsolver.options.maxtol >= settings.linsolver.options.mintol);
% Turn off undesired warnings
warning off MATLAB:bicgstab:tooSmallTolerance;
warning off MATLAB:gmres:tooSmallTolerance;
end

function defaults = mergestruct(defaults,mods)
	labels = fieldnames(mods);
	for ind = 1:length(labels)
		l = labels{ind};
		if isstruct(mods.(l)) & ~isfield(defaults,l)
			defaults.(l) = mergestruct(defaults.(l),mods.(l));
		else
			defaults.(l) = mods.(l);
		end			
	end
end
