function [B,E] = fixedblock(H,A,direction,settings)
if strcmp(settings.mode,'generic')
	[B,E] = fixedblock_generic(H,A,direction,settings);
elseif strcmp(settings.mode,'schur')
	[B,E] = fixedblock_schur(H,A,direction,settings);
elseif strcmp(settings.mode,'multicell')
	[B,E] = fixedblock_multicell(H,A,direction,settings);
elseif strcmp(settings.mode,'twosite')
	[B,E] = fixedblock_twosite(H,A,direction,settings);
else
	error(['Unrecognized mode ' settings.mode '.']);
end
end

function [B,E] = fixedblock_generic(H,A,direction,settings)
chi = size(H,1);
[D,~,d] = size(A);
eigsolver = settings.eigsolver.handle;
eigsolver_options.isreal = settings.isreal;
eigsolver_mode = 'lm';
eigsolver_options.issym = false;

if isfield(settings.advice,'B')
	eigsolver_options.v0 = reshape(settings.advice.B,[D*D*chi,1]);
end

if direction == 'l'
	% Compute right dominant eigenvector
	fapplyTv = @(v) applyTv(v,A,H,A,'l');
	[B,E] = eigsolver(fapplyTv,D^2*chi,1,eigsolver_mode,eigsolver_options);
	B = reshape(B,[D,D,chi]);
elseif direction == 'r'
	% Compute right dominant eigenvector
	fapplyTv = @(v) applyTv(v,A,H,A,'r');
	[B,E] = eigsolver(fapplyTv,D^2*chi,1,eigsolver_mode,eigsolver_options);
	B = reshape(B,[D,D,chi]);
end
end

function [B,E] = fixedblock_schur(H,A,direction,settings)
chi = size(H,1);
[D,~,d] = size(A);
id = reshape(eye([D,D]),[D*D,1]);
% Shortcuts to solvers
tol = max(eps,settings.linsolver.options.tol);
maxit = min(D*D,settings.linsolver.options.maxit);
settings.linsolver.options.v0 = [];
linsolver = settings.linsolver.handle;
eigsolver = settings.eigsolver.handle;
eigsolver_options.isreal = settings.isreal;
eigsolver_mode = 'lm';
eigsolver_options.issym = false;

if direction == 'l'
	% Calculate left block
	B = zeros([D*D,chi]);
	% Compute right dominant eigenvector of A
	if isstruct(settings) & isfield(settings.advice,'C')
		C = settings.advice.C;
		eigsolver_options.v0 = reshape((C*C').',[D*D,1]);
	end
	fapplyTv = @(v) applyTv(v,A,1,A,'r');
	[rho,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	rho = rho/(id'*rho);
	% Start from the corner
	a = chi;
	assert(H{a,a} == 1,'Element H{%d,%d} must be 1.',a,a);
	B(:,a) = id;
	for a = (chi-1):-1:1
		% Compute new Y_a
		Y = zeros([D^2,1]);
		for b = chi:-1:(a+1)
			if ~isempty(H{b,a})
				Y = Y + applyTv(B(:,b),A,H{b,a},A,'l');
			end
		end
		% Compute new B_a
		if isempty(H{a,a})
			B(:,a) = Y;
		else
			assert(isscalar(H{a,a}),'Element H{%d,%d} is not a scalar.',a,a);
			% Load initial guess
			if isstruct(settings) & isfield(settings.advice,'B')
				settings.linsolver.options.v0 = reshape(settings.advice.B(:,:,a),[D*D,1]);
			end
			if H{a,a} == 1
				applyM = @(x) x - applyTv(x,A,1,A,'l') + (rho.'*x)*id;
				b = Y - (rho.'*Y)*id;
			elseif abs(H{a,a}) < 1
				applyM = @(x) x - H{a,a}*applyTv(x,A,1,A,'l');
				b = Y;
			else
				error('Diagonal element in H has absolute value larger than 1.');
			end
			[B(:,a),~] = linsolver(applyM,Y,settings.linsolver.options);
		end
	end
	B = reshape(B,[D,D,chi]); 
	E = real(Y.'*rho);
elseif direction == 'r'
	% Calculate right block
	B = zeros([D*D,chi]);
	% Compute left dominant eigenvector of A
	if isfield(settings.advice,'C')
		C = settings.advice.C;
		eigsolver_options.v0 = reshape(C'*C,[D*D,1]);
	end
	fapplyTv = @(v) applyTv(v,A,1,A,'l');
	[rho,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	rho = rho/(id'*rho);
	% Start from the corner
	a = 1;
	assert(H{a,a} == 1,'Element H{%d,%d} must be 1.',a,a);
	B(:,a) = id;
	for a = 2:chi
		% Compute new Y_a
		Y = zeros([D^2,1]);
		for b = 1:(a-1)
			if ~isempty(H{a,b})
				Y = Y + applyTv(B(:,b),A,H{a,b},A,'r');
			end
		end
		% Compute new B_a
		if isempty(H{a,a})
			B(:,a) = Y;
		else
			assert(isscalar(H{a,a}),'Element H{%d,%d} is not a scalar.',a,a);
			% Load initial guess
			if isstruct(settings) & isfield(settings.advice,'B')
				settings.linsolver.options.v0 = reshape(settings.advice.B(:,:,a),[D*D,1]);
			end
			if H{a,a} == 1
				applyM = @(x) x - applyTv(x,A,1,A,'r') + (rho.'*x)*id;
				b = Y - (rho.'*Y)*id;
			elseif abs(H{a,a}) < 1
				applyM = @(x) x - H{a,a}*applyTv(x,A,1,A,'r');
				b = Y;
			else
				error('Diagonal element in H has absolute value larger than 1.');
			end
			[B(:,a),~] = linsolver(applyM,Y,settings.linsolver.options);
		end
	end
	E = real(Y.'*rho);
	B = reshape(B,[D,D,chi]);
else
	error(['Unrecognized direction' direction '.']);
end
end

function [B,E] = fixedblock_twosite(H,A,direction,settings)
[D,~,d] = size(A);
id = reshape(eye(D,D),[D*D,1]);
% Shortcuts to solvers
tol = max(eps,settings.linsolver.options.tol);
maxit = min(D*D,settings.linsolver.options.maxit);
settings.linsolver.options.v0 = [];
linsolver = settings.linsolver.handle;
eigsolver = settings.eigsolver.handle;
eigsolver_options.isreal = settings.isreal;
eigsolver_mode = 'lm';
eigsolver_options.issym = false;

if direction == 'l'
	% Left block
	% Compute right dominant eigenvector of A
	if isstruct(settings) & isfield(settings.advice,'C')
		C = settings.advice.C;
		eigsolver_options.v0 = reshape((C*C').',[D*D,1]);
	end
	fapplyTv = @(v) applyTv(v,A,1,A,'r');
	[rho,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	rho = rho/(id'*rho);
	% Compute left boundary
	block = ncon({A,A},{[-1,1,-2],[1,-4,-3]});
	h = ncon({conj(block),H,block},{[5,1,2,-1],[1,2,3,4],[5,3,4,-2]});
	h = reshape((h+h')/2,[D*D,1]);
	% Solve left linear system
	if isstruct(settings) & isfield(settings.advice,'B')
		settings.linsolver.options.v0 = reshape(settings.advice.B,[D*D,1]);
	end
	applyM = @(x) x - applyTv(x,A,1,A,'l') + (rho.'*x)*id;
	b = h - real(rho.'*h)*id;
	[B,~] = linsolver(applyM,b,settings.linsolver.options);
	B = reshape(B,[D,D]);
	E = real(h.'*rho);
elseif direction == 'r'
	% Right block
	% Compute left dominant eigenvector of A
	if isstruct(settings) & isfield(settings.advice,'C')
		C = settings.advice.C;
		eigsolver_options.v0 = reshape(C'*C,[D*D,1]);
	end
	fapplyTv = @(v) applyTv(v,A,1,A,'l');
	[rho,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	rho = rho/(id'*rho);
	% Compute right boundary
	block = ncon({A,A},{[-1,1,-2],[1,-4,-3]});
	h = ncon({conj(block),H,block},{[-1,1,2,5],[1,2,3,4],[-2,3,4,5]});
	h = reshape((h+h)/2,[D*D,1]);
	% Solve right linear system
	if isstruct(settings) & isfield(settings.advice,'B')
		settings.linsolver.options.v0 = reshape(settings.advice.B,[D*D,1]);
	end
	applyM = @(x) x - applyTv(x,A,1,A,'r') + (rho.'*x)*id;
	b = h - real(rho.'*h)*id;
	[B,~] = linsolver(applyM,b,settings.linsolver.options);
	B = reshape(B,[D,D]);
	E = real(rho.'*h);
else
	error(['Unrecognized direction' direction '.']);
end
end

function [B,E] = fixedblock_multicell(H,A,direction,settings)
if ~iscell(H{1})
	[B,E] = fixedblock_multicell_generic(H,A,direction,settings);
	return
end
N = length(H);
assert(isequal(size(H),size(A)),'Size mismatch in input cell arrays.');
chi = size(H{1},1);
[D,~,d] = size(A{1});
id = reshape(eye([D,D]),[D*D,1]);
B = zeros([D*D,chi]);
% Shortcuts to solvers
tol = max(eps,settings.linsolver.options.tol);
maxit = min(D*D,settings.linsolver.options.maxit);
settings.linsolver.options.v0 = [];
linsolver = settings.linsolver.handle;
eigsolver = settings.eigsolver.handle;
eigsolver_options.isreal = settings.isreal;
eigsolver_mode = 'lm';
eigsolver_options.issym = false;

Hall = combineMPO(H);

if direction == 'l'
	% Compute right dominant eigenvector of A
	if isstruct(settings) & isfield(settings.advice,'C')
		C = settings.advice.C;
		eigsolver_options.v0 = reshape((C*C').',[D*D,1]);
	end
	fapplyTv = @(v) applyTv(v,A,1,A,'r');
	[rho,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	rho = rho/(id'*rho);
	% Start from the corner
	a = chi;
	assert(all(cellfun(@(h) h == 1,Hall{a,a})),'All {%d,%d} H must be 1.',a,a);
	B(:,a) = id;
	for a = (chi-1):-1:1
		% Compute new Y_a
		Y = zeros([D^2,1]);
		for b = chi:-1:a
			if ~isempty(Hall{b,a})
				for s = 1:size(Hall{b,a},1)
					Y = Y + applyTv(B(:,b),A,Hall{b,a}(s,:),A,'l');
				end
			end
		end
		% Compute new B_a
		Hdiag = Hall{a,a}; 
		if isempty(Hdiag)
			B(:,a) = Y;
		else
			assert(all(cellfun(@(h) isscalar(h),Hdiag)),'One or more diagonal elements in H are not scalar.');
			% Load initial guess
			if isstruct(settings) & isfield(settings.advice,'B')
				settings.linsolver.options.v0 = reshape(settings.advice.B(:,:,a),[D*D,1]);
			end
			if all(cellfun(@(h) h == 1,Hdiag))
				applyM = @(x) x - applyTv(x,A,1,A,'l') + (rho.'*x)*id;
				b = Y - (rho.'*Y)*id;
			elseif all(cellfun(@(h) h < 1,Hdiag))
				applyM = @(x) x - prod(cell2mat(Hdiag))*applyTv(x,A,1,A,'l');
				b = Y;
			else
				error('One or more element in H have absolute value larger than 1.');
			end
			[B(:,a),~] = linsolver(applyM,Y,settings.linsolver.options);
		end
	end
	B = reshape(B,[D,D,chi]);
	E = real(Y.'*rho)/N;
elseif direction == 'r'
	% Compute left dominant eigenvector of A
	if isstruct(settings) & isfield(settings.advice,'C')
		C = settings.advice.C;
		eigsolver_options.v0 = reshape(C'*C,[D*D,1]);
	end
	fapplyTv = @(v) applyTv(v,A,1,A,'l');
	[rho,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	rho = rho/(id'*rho);
	% Start from the corner
	a = 1;
	assert(all(cellfun(@(h) h == 1,Hall{a,a})),'All {%d,%d} H must be 1.',a,a);
	B(:,a) = id;
	for a = 2:chi
		% Compute new Y_a
		Y = zeros([D^2,1]);
		for b = 1:(a-1)
			if ~isempty(Hall{a,b})
				for s = 1:size(Hall{a,b},1)
					Y = Y + applyTv(B(:,b),A,Hall{a,b}(s,:),A,'r');
				end
			end
		end
		% Compute new B_a
		Hdiag = Hall{a,a}; 
		if isempty(Hdiag)
			B(:,a) = Y;
		else
			assert(all(cellfun(@(h) isscalar(h),Hdiag)),'One or more diagonal elements in H are not scalar.');
			% Load initial guess
			if isstruct(settings) & isfield(settings.advice,'B')
				settings.linsolver.options.v0 = reshape(settings.advice.B(:,:,a),[D*D,1]);
			end
			if all(cellfun(@(h) h == 1,Hdiag))
				applyM = @(x) x - applyTv(x,A,1,A,'r') + (rho.'*x)*id;
				b = Y - (rho.'*Y)*id;
			elseif all(cellfun(@(h) abs(h) < 1,Hdiag))
				applyM = @(x) x - prod(cell2mat(Hdiag))*applyTv(x,A,1,A,'r');
				b = Y;
			else
				error('One or more element in H have absolute value larger than 1.');
			end
			[B(:,a),~] = linsolver(applyM,Y,settings.linsolver.options);
		end
	end
	B = reshape(B,[D,D,chi]);
	E = real(Y.'*rho)/N;
else
	error(['Unrecognized direction' direction '.']);
end
end


function [B,E] = fixedblock_multicell_generic(H,A,direction,settings)
H
N = length(H);
chi = size(H{1},1);
[D,~,d] = size(A{1});
id = reshape(eye([D,D]),[D*D,1]);
B = zeros([D*D,chi]);
% Shortcuts to solvers
tol = max(eps,settings.linsolver.options.tol);
maxit = min(D*D,settings.linsolver.options.maxit);
settings.linsolver.options.v0 = [];
linsolver = settings.linsolver.handle;
eigsolver = settings.eigsolver.handle;
eigsolver_options.isreal = settings.isreal;
eigsolver_mode = 'lm';
eigsolver_options.issym = false;

if direction == 'l'
	if isstruct(settings) & isfield(settings.advice,'H_left')
		eigsolver_options.v0 = reshape(settings.advice.H_left,[D*D*chi,1]);
	end
	% Compute right dominant eigenvector
	fapplyTv = @(v) applyTvNtimes(v,A,H,A,'l');
	[B,E] = eigsolver(fapplyTv,D^2*chi,1,eigsolver_mode,eigsolver_options);
	E = E/N;
	B = reshape(B,[D,D,chi]);
elseif direction == 'r'
	if isstruct(settings) & isfield(settings.advice,'H_right')
			eigsolver_options.v0 = reshape(settings.advice.H_right,[D*D*chi,1]);
	end
	% Compute right dominant eigenvector
	fapplyTv = @(v) applyTvNtimes(v,A,H,A,'r');
	[B,E] = eigsolver(fapplyTv,D^2*chi,1,eigsolver_mode,eigsolver_options);
	E = E/N;
	B = reshape(B,[D,D,chi]);
else
	error(['Unrecognized direction' direction '.']);
end
end
function v = applyTNtimes(M,A1,H,A2,direction)
N = length(A1);
if direction == 'l'
	ind = 1:N;
elseif direction == 'r'
	ind = N:(-1):1;
end
for n = ind
	v = applyT(M,A1{n},H{n},A2{n},direction);
end
end

function v = applyTvNtimes(v,A1,H,A2,direction)
N = length(A1);
if direction == 'l'
	ind = 1:N;
elseif direction == 'r'
	ind = N:(-1):1;
end
for n = ind
	v = applyTv(v,A1{n},H{n},A2{n},direction);
end
end
 
