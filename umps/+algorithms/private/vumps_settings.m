function settings = vumps_settings(custom_settings)
% VUMPS_SETTINGS Construct settings structure for vumps problem.

% Define default settings
settings = struct;
settings.maxit = 100;
settings.tol = 1e-12;
settings.isreal = false;
settings.verbose = true;
settings.initial = struct;
settings.advice = struct;
settings.eigsolver = struct;
settings.eigsolver.handle = @linalg.eigs;
settings.eigsolver.mode = 'sr';
settings.eigsolver.options = struct;
settings.eigsolver.options.issym = true;
settings.eigsolver.options.isreal = false;
settings.eigsolver.options.tol = eps;
settings.eigsolver.options.dynamictol = true;
settings.eigsolver.options.dynamicfactor = 1e-2;
settings.eigsolver.options.mintol = eps;
settings.eigsolver.options.maxtol = 1e-3;
settings.linsolver.handle = @linalg.bicgstab;
settings.linsolver.options = struct;
settings.linsolver.options.tol = 1e-12;
settings.linsolver.options.maxit = 100;
settings.linsolver.options.v0 = [];
settings.linsolver.options.dynamictol = true;
settings.linsolver.options.dynamicfactor = 1e-2;
settings.linsolver.options.mintol = eps;
settings.linsolver.options.maxtol = 1e-3;
% Merge settings with custom ones
if exist('custom_settings','var') && ~isempty(custom_settings)
	settings = mergestruct(settings,custom_settings);
end
% Update eigsolver settings for real matrices
if settings.isreal
	settings.eigsolver.options.isreal = true;
end

% Do some minimal checks
assert(mod(settings.maxit,1) == 0,'%s must be an integer.','settings.maxit');
assert(settings.tol <= 1 && settings.tol >= 0,'%s must be between 0 and 1.','settings.tol');
assert(settings.eigsolver.options.dynamicfactor <= 1 && settings.eigsolver.options.dynamicfactor >= 0,'%s must be between 0 and 1.','settings.eigsolver.options.dynamicfactor');
assert(settings.eigsolver.options.maxtol >= settings.eigsolver.options.mintol,'In %s, mintol must be smaller than maxtol.','settings.eigsolver.options');
assert(mod(settings.linsolver.options.maxit,1) == 0,'%s must be an integer.','settings.linsolver.options.maxit');
assert(settings.linsolver.options.tol <= 1 && settings.linsolver.options.tol >= 0,'%s must be between 0 and 1.','settings.linsolver.options.tol');
assert(settings.linsolver.options.dynamicfactor <= 1 && settings.linsolver.options.dynamicfactor >= 0,'%s must be between 0 and 1.','settings.linsolver.options.dynamicfactor');
assert(settings.linsolver.options.maxtol >= settings.linsolver.options.mintol,'In %s, mintol must be smaller than maxtol.','settings.linsolver.options');
% Turn off undesired warnings
warning off MATLAB:bicgstab:tooSmallTolerance;
warning off MATLAB:gmres:tooSmallTolerance;
end

function defaults = mergestruct(defaults,mods)
	labels = fieldnames(mods);
	for ind = 1:length(labels)
		l = labels{ind};
		if isstruct(mods.(l))
			defaults.(l) = mergestruct(defaults.(l),mods.(l));
		else
			defaults.(l) = mods.(l);
		end
	end
end
