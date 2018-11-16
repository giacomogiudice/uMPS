function [A_left,A_right,C,output,stats] = vumps_multicell(H,D,d,settings)
N = length(H);
if N == 1
	warning('VUMPS:vumps_multicell:singlecell','multicell version is suboptimal for a single site unit cell.');
end
pbc = @(n) mod(n+N-1,N) + 1;
shift = @(C,n) circshift(C,-[0,n]);

if all(isfield(settings.initial,{'A_left','A_right','C'}))
	% The initial conditions are provided
	A_left = settings.initial.A_left;
	A_right = settings.initial.A_right;
	C = settings.initial.C;
	assert(isequal(size(A_left),size(A_right),size(C)),'Size mismatch between input cell arrays.');
	assert(all(cellfun(@(Al,Ar) isequal(size(Al),size(Ar)),A_left,A_right)),'Size mismatch between left and right canonical forms.');
	assert(all(cellfun(@(Al,Ac) isequal([size(Al,1),size(Al,1)],size(Ac)),A_left,C)),'Size mismatch between canonical forms and central tensor');
	A = cell(1,N);
	% Compute error to update tolerances
	err = zeros(1,N);
	for n = 1:N
		err(n) = error_gauge([],A_left{n},A_right{n},C{pbc(n-1)},C{n});
	end
	% Update tolerances
	if settings.eigsolver.options.dynamictol
		settings.eigsolver.options.tol = update_tol(min(err),settings.eigsolver.options);
	end
	if settings.linsolver.options.dynamictol
		settings.linsolver.options.tol = update_tol(min(err),settings.linsolver.options);
	end 
	if D > size(A_left{1},1)
		% Increase the bond dimension
		B_left = cell(1,N);
		B_right = cell(1,N);
		B_left{1} = fixedblock(H,A_left,'l',settings);
		B_right{1} = fixedblock(shift(H,1),shift(A_right,1),'r',settings);
		for n = 1:(N-1)
			m = pbc(-n+2);
			B_left{n+1} = applyT(B_left{n},A_left{n},H{n},A_left{n},'l');
			B_right{pbc(m-1)} = applyT(B_right{m},A_right{m},H{m},A_right{m},'r');
		end
		for n = 1:N
			fapplyHC = @(M) applyHC(M,[],B_left{pbc(n+1)},B_right{pbc(n-1)},[],[],settings.mode);
			[Anew_left{n},Anew_right{n},Cnew{n}] = increasebond(D,A_left{n},A_right{pbc(n+1)},C{n},H{n},B_left{pbc(n+1)},B_right{pbc(n-1)});
		end
		A_right = Anew_right;
		A_left = Anew_left;
		C = Cnew;
	elseif D < size(A_left,1)
		error('Bond dimension provided is smaller than initial conditions.');
	end
	for n = 1:N
		A{n} = ncon({A_left{n},C{n}},{[-1,1,-3],[1,-2]});
	end
else 
	% Nothing is provided, generate at random
	A_left = cell(1,N);
	A_right = cell(1,N);
	A = cellfun(@(h) randn([D,D,d]) + 1i*randn([D,D,d]),H,'UniformOutput',false);
	C = cellfun(@(h) diag(rand(1,D)),H,'UniformOutput',false);
	for n = 1:N
		[A_left{n},A_right{n}] = update_canonical(A{n},C{pbc(n-1)},C{n});
	end
end
% Define solvers
eigsolver = settings.eigsolver.handle;
% Initialize stats log
stats = struct;
if nargout == 5
	savestats = true;
	stats.err = zeros(1,settings.maxit);
	stats.energy = zeros(1,settings.maxit);
	stats.energydiff = zeros(1,settings.maxit);
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
[B_left{1},energy_left] = fixedblock(H,A_left,'l',settings);
[B_right{1},energy_right] = fixedblock(shift(H,1),shift(A_right,1),'r',settings);
for n = 1:(N-1)
	m = pbc(-n+2);
	B_left{n+1} = applyT(B_left{n},A_left{n},H{n},A_left{n},'l');
	B_right{pbc(m-1)} = applyT(B_right{m},A_right{m},H{m},A_right{m},'r');
end
energy = zeros(1,N);
energy_prev = mean([energy_left,energy_right])*ones(1,N);
% Main VUMPS loop
output.flag = 1;
if settings.verbose
	fprintf('Iter\t      Energy\t Energy Diff\t Gauge Error\tLap Time [s]\n')
	fprintf('   0\t%12g\n',mean(energy_prev));
end
for iter = 1:settings.maxit
	tic
	for n = 1:N
		% Solve effective problem for A
		settings.eigsolver.options.v0 = reshape(A{n},[D*D*d,1]);
		applyHAv = @(v) reshape(applyHA(reshape(v,[D,D,d]),H{n},B_left{n},B_right{n},[],[],settings.mode),[D*D*d,1]);
		[Av,~] = eigsolver(applyHAv,D*D*d,1,settings.eigsolver.mode,settings.eigsolver.options);
		A{n} = reshape(Av,[D,D,d]);
		% Solve effective problem for C_left
		B_mid = applyT(B_right{n},A_right{n},H{n},A_right{n},'r');
		settings.eigsolver.options.v0 = reshape(C{pbc(n-1)},[D*D,1]);
		applyHCv = @(v) reshape(applyHC(reshape(v,[D,D]),[],B_left{n},B_mid,[],[],settings.mode),[D*D,1]);
		[Cv,~] = eigsolver(applyHCv,D*D,1,settings.eigsolver.mode,settings.eigsolver.options);
		Cv = Cv/sign(Cv(1));
		C_left = reshape(Cv,[D,D]);
		% Solve effective problem for C_right
		B_mid = applyT(B_left{n},A_left{n},H{n},A_left{n},'l');
		settings.eigsolver.options.v0 = reshape(C{n},[D*D,1]);
		applyHCv = @(v) reshape(applyHC(reshape(v,[D,D]),[],B_mid,B_right{n},[],[],settings.mode),[D*D,1]);
		[Cv,~] = eigsolver(applyHCv,D*D,1,settings.eigsolver.mode,settings.eigsolver.options);
		Cv = Cv/sign(Cv(1));
		C_right = reshape(Cv,[D,D]);
		% Update the canonical forms
		[A_left{n},A_right{n}] = update_canonical(A{n},C_left,C_right);
		C{pbc(n-1)} = C_left;
		C{n} = C_right;
		% Get error
		err(n) = error_gauge(A{n},A_left{n},A_right{n},C_left,C_right);
		% Update the next environment blocks
		settings.advice.C = C{n};
		settings.advice.B = B_left{pbc(n+1)};
		[B_left{pbc(n+1)},energy_left] = fixedblock(shift(H,n),shift(A_left,n),'l',settings);
		settings.advice.C = C{pbc(n+1)};
		settings.advice.B = B_right{pbc(n+1)};
		[B_right{pbc(n+1)},energy_right] = fixedblock(shift(H,n+1),shift(A_right,n+1),'r',settings);
		energy(n) = mean([energy_left,energy_right]);
	end
	laptime = toc;
	% Update tolerances
	if settings.eigsolver.options.dynamictol
		settings.eigsolver.options.tol = update_tol(min(err),settings.eigsolver.options);
	end
	if settings.linsolver.options.dynamictol
		settings.linsolver.options.tol = update_tol(min(err),settings.linsolver.options);
	end
	% Print results of interation
	if settings.verbose
		fprintf('%4d\t%12g\t%12g\t%12g%12.1f\n',iter,mean(energy),mean(energy_prev - energy),max(err),laptime);
	end
	if savestats
		stats.err(iter) = max(err);
		stats.energy(iter) = mean(energy);
		stats.energydiff(iter) = mean(energy_prev - energy);
	end
	if max(err) < settings.tol
		output.flag = 0;
		break
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
end
% Choose gauge in which each C is diagonal
for n = 1:N
	[U,S,V] = svd(C{n},'econ');
	A_left{n} = ncon({A_left{n},U},{[-1,2,-3],[2,-2]});
	A_left{pbc(n+1)} = ncon({U',A_left{pbc(n+1)}},{[-1,1],[1,-2,-3]});
	A_right{n} = ncon({A_right{n},V},{[-1,1,-3],[1,-2]});
	A_right{pbc(n+1)} = ncon({V',A_right{pbc(n+1)}},{[-1,1],[1,-2,-3]});
	C{n} = S;
end
end

