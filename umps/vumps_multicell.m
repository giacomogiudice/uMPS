function [A_left,A_right,C,output,stats] = vumps_multicell(H,D,d,settings)
N = length(H);
fixedblock = @fixedblock_m;
applyHA = @applyHA_s;
applyHC = @applyHC_g;
applyH2s = @applyH2s_s;
pbc = @(n) mod(n+N-1,N) + 1;
% Intialize at random
A_left = cell(1,N);
A_right = cell(1,N);
C = cell(1,N);;

for n = 1:N
	A = randn([D,D,d]) + 1i*randn([D,D,d]);
	C{n} = diag(rand(1,D));
	[A_left{n},A_right{n}] = update_canonical(A,C{n});
end

% Define solvers
eigsolver = settings.eigsolver.handle;
% settings.eigsolver.options.v0 = [];
% Initialize stats log
stats = struct;
if nargout == 5
	savestats = true;
	stats.err = zeros(1,settings.maxit);
	stats.energy = zeros(1,settings.maxit);
	stats.energydiff = zeros(1,settings.maxit);
end
% TODO Update tolerances
% Main VUMPS loop
output.flag = 1;
if settings.verbose
	fprintf('Iter\t      Energy\t Energy Diff\t Gauge Error\tLap Time [s]\n')
	% fprintf('   0\t%12g\n',energy_prev);
end
energy_prev = inf;
for iter = 1:settings.maxit
	tic
	for n = 1:N
		% Update environments
		[L,energy_left] = fixedblock(circshift(H,-[0,n-1]),circshift(A_left,-[0,n-1]),circshift(A_right,-[0,n-1]),'l',settings);
		[R,energy_right] = fixedblock(circshift(H,-[0,n]),circshift(A_left,-[0,n]),circshift(A_right,-[0,n]),'r',settings);
		L_prime = applyT(L,A_left{n},H{n},A_left{n},'l');
		R_prime = applyT(R,A_right{n},H{n},A_right{n},'r');
		% Solve effective problem for A
		% settings.eigsolver.options.v0 = reshape(A,[D*D*d,1]);
		applyHAv = @(v) reshape(applyHA(reshape(v,[D,D,d]),H{n},L,R),[D*D*d,1]);
		[Av,~] = eigsolver(applyHAv,D*D*d,1,settings.eigsolver.mode,settings.eigsolver.options);
		A = reshape(Av,[D,D,d]);
		% Solve effective problem for C_left
		% settings.eigsolver.options.v0 = reshape(C{n},[D*D,1]);
		applyHCv = @(v) reshape(applyHC(reshape(v,[D,D]),[],L,R_prime),[D*D,1]);
		[Cv,~] = eigsolver(applyHCv,D*D,1,settings.eigsolver.mode,settings.eigsolver.options);
		C_left = reshape(Cv,[D,D]);
		% Solve effective problem for C_right
		% settings.eigsolver.options.v0 = reshape(C{pbc(n)},[D*D,1]);
		applyHCv = @(v) reshape(applyHC(reshape(v,[D,D]),[],L_prime,R),[D*D,1]);
		[Cv,~] = eigsolver(applyHCv,D*D,1,settings.eigsolver.mode,settings.eigsolver.options);
		C_right = reshape(Cv,[D,D]);
		% Update the canonical forms
		[A_left{n},A_right{n}] = update_canonical_m(A,C_left,C_right);
		C{pbc(n-1)} = C_right;
		C{n} = C_left;
		% Get error
		err(n) = error_gauge_m(A,C_left,C_right,A_left{n},A_right{n});
		energy(n) = mean([energy_left,energy_right]);
		% Update the environment blocks
		laptime = toc;
	end
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
% Set output information
output.energy = mean(energy);
output.iter = iter;
output.err = max(err);
if savestats
	stats.err = stats.err(1:iter);
	stats.energy = stats.energy(1:iter);
	stats.energydiff = stats.energydiff(1:iter);
end
% Choose gauge in which each C is diagonal
for n = 1:N
	[U,S,V] = svd(C{n},'econ');
	A_left{n} = ncon({A_left{n},U},{[-1,2,-3],[2,-2]});
	A_right{pbc(n+1)} = ncon({V',A_right{pbc(n+1)}},{[-1,1],[1,-2,-3]});
	C{n} = S;
end