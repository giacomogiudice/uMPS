function [A_left,A_right,C,A,output,blocks,stats] = vumps_multicell(H,D_list,d,settings)
N = length(H);
if N == 1
	warning('VUMPS:vumps_multicell:singlecell','multicell version is suboptimal for a single site unit cell.');
end
% Function to wrap around the indices
pbc = @(n) mod(n+N-1,N) + 1;
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
	assert(isequal(size(A_left),size(A_right),size(C),size(H)),'Size mismatch between input cell arrays.');
	assert(all(cellfun(@(Al,Ar) isequal(size(Al),size(Ar)),A_left,A_right)),'Size mismatch between left and right canonical forms.');
	assert(all(cellfun(@(Al,Ac) isequal([size(Al,1),size(Al,1)],size(Ac)),A_left,C)),'Size mismatch between canonical forms and central tensor');
	A = cell(1,N);
	for n = 1:N
		A{n} = ncon({A_left{n},C{n}},{[-1,1,-3],[1,-2]});
	end
	D = size(C{1},1);
	assert(D <= D_list(1),'Bond dimension of initial MPS must be smaller or equal to first element of D_list.');
else 
	% Nothing is provided, generate at random
	A_left = cell(1,N);
	A_right = cell(1,N);
	A = cell(1,N);
	A = cellfun(@(a) randn([D,D,d]),A,'UniformOutput',false);
	if ~settings.isreal
		A = cellfun(@(a) a + 1i*randn([D,D,d]),A,'UniformOutput',false);
	end
	C = cellfun(@(h) diag(rand(1,D)),H,'UniformOutput',false);
	for n = 1:N
		[A_left{n},A_right{n}] = update_canonical(A{n},C{pbc(n-1)},C{n});
	end
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
for n = 1:N
	err(n) = error_gauge(A{n},A_left{n},A_right{n},C{pbc(n-1)},C{n});
end
% Update tolerances
if settings.eigsolver.options.dynamictol
	settings.eigsolver.options.tol = update_tol(min(err),settings.eigsolver.options);
end
if settings.linsolver.options.dynamictol
	settings.linsolver.options.tol = update_tol(min(err),settings.linsolver.options);
end
% Generate the environment blocks
B_left = cell(1,N);
B_right = cell(1,N);
[B_left{1},B_right{1},energy] = update_environments(A_left,C{1},A_right,C{2},H,B_left,B_right,settings);
for n = 1:(N-1)
	m = pbc(-n+2);
	B_left{n+1} = applyT(B_left{n},A_left{n},H{n},A_left{n},'l');
	B_right{pbc(m-1)} = applyT(B_right{m},A_right{m},H{m},A_right{m},'r');
end
energy_prev = energy*ones(1,N);
energy = zeros(1,N);
% Main VUMPS loop
output.flag = 2;
if settings.verbose
	fprintf('Iter\t      Energy\t Energy Diff\t Gauge Error\t    Bond Dim\tLap Time [s]\n')
	fprintf('   0\t%12g\n',mean(energy_prev));
end
for iter = 1:settings.maxit
	tic
	for n = 1:N
		% Solve effective problems for A, C_left and C_right
		[A{n},C_left,C_right] = ...
			solve_local(A_left{n},C{pbc(n-1)},A_right{n},C{n},A{n},H{n},B_left{n},B_right{n},settings);
		% Update the canonical forms
		[A_left{n},A_right{n}] = update_canonical(A{n},C_left,C_right);
		C{pbc(n-1)} = C_left;
		C{n} = C_right;
		% Update the next environment blocks
		[B_left{pbc(n+1)},B_right{pbc(n+1)},energy(n)] = ...
			update_environments(shift(A_left,n),C{n},shift(A_right,n),C{pbc(n+1)},shift(H,n),B_left{pbc(n+1)},B_right{pbc(n+1)},settings);
		% Get error
		err(n) = error_gauge(A{n},A_left{n},A_right{n},C_left,C_right);
	end
	laptime = toc;
	% Update tolerances
	if settings.eigsolver.options.dynamictol
		settings.eigsolver.options.tol = update_tol(min(err),settings.eigsolver.options);
	end
	if settings.linsolver.options.dynamictol
		settings.linsolver.options.tol = update_tol(min(err),settings.linsolver.options);
	end
	energydiff = energy - energy_prev;
	% Print results of interation
	if settings.verbose
		fprintf('%4d\t%12g\t%12g\t%12g%12d%12.1f\n',iter,mean(energy),mean(energydiff),max(err),D_list(bond_ind),laptime);
	end
	if savestats
		stats.err(iter) = max(err);
		stats.energy(iter) = mean(energy);
		stats.energydiff(iter) = mean(energydiff);
		stats.bond(iter) = D;
	end
	if bond_ind == length(D_list)
		if max(err) < settings.tol
			% Stopping condition on the gauge error
			output.flag = 0;
			break
		elseif mean(abs(energy_prev - energy)) < eps
			% Stopping condition on stagnation
			output.flag = 1;
			break
		end
	elseif err < growtol(bond_ind)
		% Increase bond dimension
		bond_ind = bond_ind + 1;
		D = D_list(bond_ind);
		for n = 1:N
			[Anew_left{n},Anew_right{pbc(n+1)},Cnew{n},Anew{n},B_left{pbc(n+1)},B_right{pbc(n-1)}] = ...
				increasebond(D,A_left{n},A_right{pbc(n+1)},C{n},H{n},B_left{pbc(n+1)},B_right{pbc(n-1)});
		end
		A_right = Anew_right;
		A_left = Anew_left;
		C = Cnew;
		A = Anew;
		[B_left{1},B_right{1}] = update_environments(A_left,C{1},A_right,C{2},H,B_left,B_right,settings);
		for n = 1:(N-1)
			m = pbc(-n+2);
			B_left{n+1} = applyT(B_left{n},A_left{n},H{n},A_left{n},'l');
			B_right{pbc(m-1)} = applyT(B_right{m},A_right{m},H{m},A_right{m},'r');
		end
	end
	energy_prev = energy;
end
% Update blocks
for n = 1:(N-1)
	m = pbc(-n+2);
	B_left{n+1} = applyT(B_left{n},A_left{n},H{n},A_left{n},'l');
	B_right{pbc(m-1)} = applyT(B_right{m},A_right{m},H{m},A_right{m},'r');
end
% Compute estimate of variance
g = zeros(1,N);
for n = 1:N
	g(n) = error_variance(A_left{n},C{n},A_right{pbc(n+1)},{H{n},H{pbc(n+1)}},B_left{n},B_right{pbc(n+1)});
end
% Set output information
output.energy = mean(energy);
output.iter = iter;
output.err = max(err);
output.energyvariance = mean(g);

if savestats
	stats.err = stats.err(1:iter);
	stats.energy = stats.energy(1:iter);
	stats.energydiff = stats.energydiff(1:iter);
	stats.bond = stats.bond(1:iter);
end
% Choose gauge in which each C is diagonal
for n = 1:N
	[U,S,V] = svd(C{n},'econ');
	A_left{n} = ncon({A_left{n},U},{[-1,2,-3],[2,-2]});
	A_left{pbc(n+1)} = ncon({U',A_left{pbc(n+1)}},{[-1,1],[1,-2,-3]});
	A_right{n} = ncon({A_right{n},V},{[-1,1,-3],[1,-2]});
	A_right{pbc(n+1)} = ncon({V',A_right{pbc(n+1)}},{[-1,1],[1,-2,-3]});
	A{n} = ncon({A{n},V},{[-1,1,-3],[1,-2]});
	A{pbc(n+1)} = ncon({U',A{pbc(n+1)}},{[-1,1],[1,-2,-3]});
	C{n} = S;
end
% Store other extra
if nargout >= 5
	[B_left{1},B_right{1}] = update_environments(A_left,C{1},A_right,C{2},H,B_left,B_right,settings);
	for n = 1:(N-1)
		m = pbc(-n+2);
		B_left{n+1} = applyT(B_left{n},A_left{n},H{n},A_left{n},'l');
		B_right{pbc(m-1)} = applyT(B_right{m},A_right{m},H{m},A_right{m},'r');
	end
	blocks.left = B_left;
	blocks.right = B_right;
end
end

function D = shift(C,n)
D = circshift(C,-[0,n]);
end

function [A,C_left,C_right] = solve_local(A_left,C_left,A_right,C_right,A,H,B_left,B_right,settings)
% Notice that this function takes the tensors themselves and not the cell arrays
[D,~,d] = size(A);
% Define solvers
eigsolver = settings.eigsolver.handle;
if settings.isreal && isequal(settings.eigsolver.mode,'sr')
	settings.eigsolver.mode = 'sa';
end
% Solve effective problem for A
settings.eigsolver.options.v0 = reshape(A,[D*D*d,1]);
applyHAv = @(v) reshape(applyHA(reshape(v,[D,D,d]),H,B_left,B_right,[],[],settings.mode),[D*D*d,1]);
[Av,~] = eigsolver(applyHAv,D*D*d,1,settings.eigsolver.mode,settings.eigsolver.options);
A = reshape(Av,[D,D,d]);
% Solve effective problem for C_left
B_mid = applyT(B_right,A_right,H,A_right,'r');
settings.eigsolver.options.v0 = reshape(C_left,[D*D,1]);
applyHCv = @(v) reshape(applyHC(reshape(v,[D,D]),[],B_left,B_mid,[],[],settings.mode),[D*D,1]);
[Cv,~] = eigsolver(applyHCv,D*D,1,settings.eigsolver.mode,settings.eigsolver.options);
Cv = Cv/sign(Cv(1));
C_left = reshape(Cv,[D,D]);
% Solve effective problem for C_right
B_mid = applyT(B_left,A_left,H,A_left,'l');
settings.eigsolver.options.v0 = reshape(C_right,[D*D,1]);
applyHCv = @(v) reshape(applyHC(reshape(v,[D,D]),[],B_mid,B_right,[],[],settings.mode),[D*D,1]);
[Cv,~] = eigsolver(applyHCv,D*D,1,settings.eigsolver.mode,settings.eigsolver.options);
Cv = Cv/sign(Cv(1));
C_right = reshape(Cv,[D,D]);
end

function [B_left,B_right,energy] = update_environments(A_left,C_left,A_right,C_right,H,B_left,B_right,settings)
if ~isempty(C_left)
	settings.advice.C = C_left;
end
if ~isempty(B_left) && size(B_left,1) == size(C_left,1)
	settings.advice.B = B_left;
end
[B_left,energy_left] = fixedblock(H,A_left,'l',settings);
if ~isempty(C_right)
	settings.advice.C = C_right;
end
if ~isempty(B_right) && size(B_right,1) == size(C_right,1)
	settings.advice.B = B_right;
end
[B_right,energy_right] = fixedblock(shift(H,1),shift(A_right,1),'r',settings);
energy = mean([energy_left,energy_right]);
% Fix normalization in case of generic MPO
if strcmp(settings.mode,'generic')
	B_norm = sqrt(abs(ncon({B_left,conj(C_left),C_left,B_right},{[1,4,3],[1,2],[4,5],[2,5,3]})));
	B_left = B_left/B_norm;
	B_right = B_right/B_norm;
	energy = real(energy);
end
end

