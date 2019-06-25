function [A_left,A_right,C,A,output,blocks,stats] = vumps(H,D_list,d,settings)
% Parse settings
if ~exist('settings','var') || isempty(settings)
	settings = vumps_settings();
else
	settings = vumps_settings(settings);
end
% Select mode of operation
if isequal(settings.mode,'schur')
	chi = size(H,1);
	assert(iscell(H),'Input operator should be a cell array.');
	assert(isequal(size(H),[chi,chi]),'Size mismatch in input operator.');
elseif isequal(settings.mode,'generic')
	chi = size(H,1);
	assert(isequal(size(H),[chi,chi,d,d]),'Input operator should be a rank-4 tensor.');
elseif isequal(settings.mode,'twosite')
	assert(isequal(size(H),[d,d,d,d]),'Input operator should be a rank-4 tensor.');
elseif isequal(settings.mode,'multicell')
	 [A_left,A_right,C,A,output,blocks,stats] = vumps_multicell(H,D_list,d,settings);
	return
else
	error(['Unrecognized mode ' settings.mode '.'])
end

% Do checks on D_list
assert(all(diff(D_list) > 0),'D_list must be a vector of positive integers in ascending order.');
D = D_list(1);
bond_ind = 1;
growtol = logspace(0,log10(settings.tol),length(D_list)+1);
growtol(1) = [];
growtol(end) = 0;

if all(isfield(settings.initial,{'A_left','A_right','C'}))
	% The initial conditions are provided
	A_left = settings.initial.A_left;
	A_right = settings.initial.A_right;
	C = settings.initial.C;
	assert(isequal(size(A_left),size(A_right)),'Size mismatch between left and right canonical forms.');
	assert(isequal([size(A_left,1),size(A_left,2)],size(C)),'Size mismatch between canonical forms and central tensor.');
	A = ncon({A_left,C},{[-1,1,-3],[1,-2]});
	D = size(C,1);
	assert(D <= D_list(1),'Bond dimension of initial MPS must be smaller or equal to first element of D_list.');
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
% Initialize stats log
stats = struct;
savestats = false;
if nargout >= 6
	savestats = true;
	stats.err = zeros(1,settings.maxit);
	stats.energy = zeros(1,settings.maxit);
	stats.energydiff = zeros(1,settings.maxit);
	stats.bond = zeros(1,settings.maxit);
end
% Build left and right blocks
err = error_gauge(A,A_left,A_right,C);
% Update tolerances
if settings.eigsolver.options.dynamictol
	settings.eigsolver.options.tol = update_tol(err,settings.eigsolver.options);
end
if settings.linsolver.options.dynamictol
	settings.linsolver.options.tol = update_tol(err,settings.linsolver.options);
end
% Generate the environment blocks
[B_left,B_right,energy_prev] = update_environments(A_left,C,A_right,H,[],[],settings);
% Main VUMPS loop
output.flag = 2;
if settings.verbose
	fprintf('Iter\t      Energy\t Energy Diff\t Gauge Error\t    Bond Dim\tLap Time [s]\n')
	fprintf('   0\t%12g\n',energy_prev);
end

for iter = 1:settings.maxit
	tic
	% Solve effective problems for A and C
	if ~all([isreal(B_left),isreal(B_right)])
		settings.isreal = 0;
		settings.eigsolver.options.isreal = 0;
	end
	[A,C] = solve_local(A_left,C,A_right,A,H,B_left,B_right,settings);
	% Update the canonical forms
	[A_left,A_right] = update_canonical(A,C);
	% Update the environment blocks
	[B_left,B_right,energy] = update_environments(A_left,C,A_right,H,B_left,B_right,settings);
	% Update tolerances
	if settings.eigsolver.options.dynamictol
		settings.eigsolver.options.tol = update_tol(err,settings.eigsolver.options);
	end
	if settings.linsolver.options.dynamictol
		settings.linsolver.options.tol = update_tol(err,settings.linsolver.options);
	end
	laptime = toc;
	% Get error
	err = error_gauge(A,A_left,A_right,C);
	% Print results of interation
	if settings.verbose
		fprintf('%4d\t%12g\t%12g\t%12g%12d%12.1f\n',iter,energy,energy_prev - energy,err,D_list(bond_ind),laptime);
	end
	if savestats
		stats.err(iter) = err;
		stats.energy(iter) = energy;
		stats.energydiff(iter) = energy_prev - energy;
		stats.bond(iter) = D;
	end	
	if bond_ind == length(D_list)
		if err < settings.tol
			% Stopping condition on the gauge error
			output.flag = 0;
			break
		elseif abs(energy_prev - energy) < eps
			% Stopping condition on stagnation
			output.flag = 1;
			break
		end
	elseif err < growtol(bond_ind)
		% Increase bond dimension
		bond_ind = bond_ind + 1;
		D = D_list(bond_ind);
		[A_left,A_right,C,A] = increasebond(D,A_left,A_right,C,H,B_left,B_right);
		[B_left,B_right] = update_environments(A_left,C,A_right,H,B_left,B_right,settings);
	end
	energy_prev = energy;
end
% Set output information
output.iter = iter;
output.err = err;
output.energy = energy;
output.energyvariance = error_variance(A_left,C,A_right,H,B_left,B_right);
if savestats
	stats.err = stats.err(1:iter);
	stats.energy = stats.energy(1:iter);
	stats.energydiff = stats.energydiff(1:iter);
	stats.bond = stats.bond(1:iter);
end
% Choose gauge in which C is diagonal
[U,S,V] = svd(C,'econ');
A_left = ncon({U',A_left,U},{[-1,1],[1,2,-3],[2,-2]});
A_right = ncon({V',A_right,V},{[-1,1],[1,2,-3],[2,-2]});
A = ncon({U',A,V},{[-1,1],[1,2,-3],[2,-2]});
C = S;

if nargout >= 5
	[B_left,B_right] = update_environments(A_left,C,A_right,H,B_left,B_right,settings);
	blocks.left =  B_left;
	blocks.right = B_right;
end
end


function [A,C] = solve_local(A_left,C,A_right,A,H,B_left,B_right,settings)
[D,~,d] = size(A);
% Define solvers
eigsolver = settings.eigsolver.handle;
if settings.isreal && isequal(settings.eigsolver.mode,'sr')
	settings.eigsolver.mode = 'sa';
end
% Solve effective problem for A
settings.eigsolver.options.v0 = reshape(A,[D*D*d,1]);
applyHAv = @(v) reshape(applyHA(reshape(v,[D,D,d]),H,B_left,B_right,A_left,A_right,settings.mode),[D*D*d,1]);
[Av,~] = eigsolver(applyHAv,D*D*d,1,settings.eigsolver.mode,settings.eigsolver.options);
A = reshape(Av,[D,D,d]);
% Solve effective problem for C
settings.eigsolver.options.v0 = reshape(C,[D*D,1]);
applyHCv = @(v) reshape(applyHC(reshape(v,[D,D]),H,B_left,B_right,A_left,A_right,settings.mode),[D*D,1]);
[Cv,~] = eigsolver(applyHCv,D*D,1,settings.eigsolver.mode,settings.eigsolver.options);
C = reshape(Cv,[D,D]);
end

function [B_left,B_right,energy] = update_environments(A_left,C,A_right,H,B_left,B_right,settings)
if ~isempty(C)
	settings.advice.C = C;
end
if ~isempty(B_left) && size(B_left,1) == size(C,1)
	settings.advice.B = B_left;
end
[B_left,energy_left] = fixedblock(H,A_left,'l',settings);
if ~isempty(B_right) && size(B_right,1) == size(C,1)
	settings.advice.B = B_right;
end
[B_right,energy_right] = fixedblock(H,A_right,'r',settings);
energy = mean([energy_left,energy_right]);
% Fix normalization in case of generic MPO
if strcmp(settings.mode,'generic')
	B_norm = sqrt(abs(ncon({B_left,B_right},{[1,2,3],[1,2,3]})));
	B_left = B_left/B_norm;
	B_right = B_right/B_norm;
	energy = real(energy);
end
end
