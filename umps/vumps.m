function [A_left,A_right,C,output,stats] = vumps(H,D,d,settings)
% Parse settings
if nargin == 3
	settings = vumps_settings();
else
	settings = vumps_settings(settings);
end
% Select mode of operation
if isequal(settings.mode,'schur')
	chi = size(H,1);
	assert(iscell(H));
	assert(isequal(size(H),[chi,chi]));
	fixedblocks = @fixedblocks_s;
	applyHA = @applyHA_s;
	applyHC = @applyHC_g;
	applyH2s = @applyH2s_s;
elseif isequal(settings.mode,'generic')
	chi = size(H,1);
	assert(isequal(size(H),[chi,chi]));
	fixedblocks = @fixedblocks_g;
	applyHA = @applyHA_g;
	applyHC = @applyHC_g;
	applyH2s = @applyH2s_g;
elseif isequal(settings.mode,'twosite')
	assert(isequal(size(H),[d,d,d,d]));
	fixedblocks = @fixedblocks_t;
	applyHA = @applyHA_t;
	applyHC = @applyHC_t;
	applyH2s = @applyH2s_t;
elseif isequal(settings.mode,'multicell')
	error('Multicell mode not available.');
else
	error(['Unrecognized mode' settings.mode])
end

if all(isfield(settings.initial,{'A_left','A_right','C'}))
	% The initial conditions are provided
	A_left = settings.initial.A_left;
	A_right = settings.initial.A_right;
	C = settings.initial.C;
	A = ncon({A_left,C},{[-1,1,-3],[1,-2]});
	assert(isequal(size(A_left),size(A_right)));
	assert(isequal([size(A_left,1),size(A_left,2)],size(C)));
	% Compute error to update tolerances
	err = max(error_gauge(A,C,A_left,A_right));
	% Update tolerances
	if settings.eigsolver.options.dynamictol
		settings.eigsolver.options.tol = update_tol(err,settings.eigsolver.options);
	end
	if settings.linsolver.options.dynamictol
		settings.linsolver.options.tol = update_tol(err,settings.linsolver.options);
	end 
	if D > size(A_left,1)
		% Increase the bond dimension
		[H_left,H_right] = fixedblocks(H,A_left,A_right,settings);
		fapplyHC = @(M) applyHC(M,H,H_left,H_right,A_left,A_right);
		[A_left,A_right,C] = increasebond(D,A_left,A_right,C,fapplyHC);
		A = ncon({A_left,C},{[-1,1,-3],[1,-2]});
	elseif D < size(A_left,1)
		error('Bond dimension provided is smaller than initial conditions');
	end
else 
	% Nothing is provided, generate at random
	A = randn([D,D,d]);
	if ~settings.isreal
		A = A + 1i*randn([D,D,d]);
	end
	C = diag(rand(D,1));
	C = C/trace(C*C');
	[A_left,A_right] = update_canonical(A,C);
end
% Define solvers
eigsolver = settings.eigsolver.handle;
% Initialize stats log
stats = struct;
savestats = false;
if nargout == 5
	savestats = true;
	stats.err = zeros(1,settings.maxit);
	stats.energy = zeros(1,settings.maxit);
	stats.energydiff = zeros(1,settings.maxit);
end
% Build left and right blocks
err = max(error_gauge(A,C,A_left,A_right));
% Update tolerances
if settings.eigsolver.options.dynamictol
	settings.eigsolver.options.tol = update_tol(err,settings.eigsolver.options);
end
if settings.linsolver.options.dynamictol
	settings.linsolver.options.tol = update_tol(err,settings.linsolver.options);
end
% Update the environment blocks
[H_left,H_right,energy_prev] = fixedblocks(H,A_left,A_right,settings);
% Main VUMPS loop
output.flag = 1;
if settings.verbose
	fprintf('Iter\t      Energy\t Energy Diff\t Gauge Error\n')
	fprintf('   0\t%12g\n',energy_prev);
end
for iter = 1:settings.maxit
	% Solve effective problem for A
	settings.eigsolver.options.v0 = reshape(A,[D*D*d,1]);
	applyHAv = @(v) reshape(applyHA(reshape(v,[D,D,d]),H,H_left,H_right,A_left,A_right),[D*D*d,1]);
	[Av,~] = eigsolver(applyHAv,D*D*d,1,settings.eigsolver.mode,settings.eigsolver.options);
	A = reshape(Av,[D,D,d]);
	% Solve effective problem for C
	settings.eigsolver.options.v0 = reshape(C,[D*D,1]);
	applyHCv = @(v) reshape(applyHC(reshape(v,[D,D]),H,H_left,H_right,A_left,A_right),[D*D,1]);
	[Cv,~] = eigsolver(applyHCv,D*D,1,settings.eigsolver.mode,settings.eigsolver.options);
	C = reshape(Cv,[D,D]);
	% Update the canonical forms
	[A_left,A_right] = update_canonical(A,C);
	% Update advice for next calculations
	settings.advice.C = C;
	settings.advice.H_left = H_left;
	settings.advice.H_right = H_right;
	% Update the environment blocks
	[H_left,H_right,energy] = fixedblocks(H,A_left,A_right,settings);
	% Get error
	err = max(error_gauge(A,C,A_left,A_right));
	if settings.verbose
		fprintf('%4d\t%12g\t%12g\t%12g\n',iter,mean(energy),mean(energy_prev - energy),max(err));
	end
	if savestats
		stats.err(iter) = err;
		stats.energy(iter) = energy;
		stats.energydiff(iter) = energy_prev - energy;
	end
	if err < settings.tol
		output.flag = 0;
		break
	end
	% Update tolerances
	if settings.eigsolver.options.dynamictol
		settings.eigsolver.options.tol = update_tol(err,settings.eigsolver.options);
	end
	if settings.linsolver.options.dynamictol
		settings.linsolver.options.tol = update_tol(err,settings.linsolver.options);
	end
	energy_prev = energy;
end
% Set output information
output.iter = iter;
output.err = err;
output.energy = energy;
fapplyHA = @(M) applyHA(M,H,H_left,H_right,A_left,A_right);
fapplyH2s = @(M) applyH2s(M,H,H_left,H_right);
output.energyvariance = error_variance(A_left,A_right,C,fapplyHA,fapplyH2s);
if savestats
	stats.err = stats.err(1:iter);
	stats.energy = stats.energy(1:iter);
	stats.energydiff = stats.energydiff(1:iter);
end
% Choose gauge in which C is diagonal
[U,S,V] = svd(C,'econ');
A_left = ncon({U',A_left,U},{[-1,1],[1,2,-3],[2,-2]});
A_right = ncon({V',A_right,V},{[-1,1],[1,2,-3],[2,-2]});
C = S;
end