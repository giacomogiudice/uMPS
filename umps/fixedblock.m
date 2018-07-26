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

if direction == 'l'
	% Compute right dominant eigenvector
	fapplyTv = @(v) applyTv(v,A,H,A,'l');
	[B,E] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	B = reshape(B,[D,D,chi]);
elseif direction == 'r'
	% Compute right dominant eigenvector
	fapplyTv = @(v) applyTv(v,A,H,A,'r');
	[B,E] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
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
	if H{a,a} ~= 1
		error(['Element H{' num2str(a) ',' num2str(a) '} must be 1.']);
	end
	B(:,a) = id;
	for a = (chi-1):-1:1
		% Compute new Y_a
		Y = zeros([D^2,1]);
		for b = chi:-1:(a+1)
			if ~isempty(H{b,a})
				Y = Y + applyTv(B(:,b),A,H{b,a},A,'l');
			end
		end
		% Load initial guess
		if isfield(settings.advice,'B')
			settings.linsolver.options.v0 = reshape(settings.advice.B(:,:,a),[D*D,1]);
		end
		% Compute new L_a
		if ~isempty(H{a,a}) & ~isscalar(H{a,a})
			error(['Element H{' num2str(a) ',' num2str(a) '} is not a scalar.']);
		elseif isempty(H{a,a}) || H{a,a} == 0
			B(:,a) = Y;
		elseif H{a,a} == 1
			applyM = @(x) x - applyTv(x,A,1,A,'l') + (rho.'*x)*id;
			[B(:,a),~] = linsolver(applyM,Y - (rho.'*Y)*id,settings.linsolver.options);
		elseif abs(H{a,a}) < 1
			applyM = @(x) x - H{a,a}*applyTv(x,A,1,A,'l');
			[B(:,a),~] = linsolver(applyM,Y,settings.linsolver.options);
		else
			error(['Element H{' num2str(a) ',' num2str(a) '} has absolute value larger than 1.']);
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
	if H{a,a} ~= 1
		error(['Element H{' num2str(a) ',' num2str(a) '} must be 1.']);
	end
	B(:,a) = id;
	for a = 2:chi
		% Compute new Y_a
		Y = zeros([D^2,1]);
		for b = 1:(a-1)
			if ~isempty(H{a,b})
				Y = Y + applyTv(B(:,b),A,H{a,b},A,'r');
			end
		end
		% Load initial guess
		if isstruct(settings) & isfield(settings.advice,'B')
			settings.linsolver.options.v0 = reshape(settings.advice.B(:,:,a),[D*D,1]);
		end
		% Compute new R_a 
		if ~isempty(H{a,a}) & ~isscalar(H{a,a})
			error(['Element H{' num2str(a) ',' num2str(a) '} is not a scalar.']);
		elseif isempty(H{a,a}) || H{a,a} == 0
			B(:,a) = Y;
		elseif H{a,a} == 1
			applyM = @(x) x - applyTv(x,A,1,A,'r') + (rho.'*x)*id;
			[B(:,a)] = linsolver(applyM,Y - (rho.'*Y)*id,settings.linsolver.options);
		elseif abs(H{a,a}) < 1
			applyM = @(x) x - H{a,a}*applyTv(x,A,1,A,'r');
			[B(:,a),~] = linsolver(applyM,Y,settings.linsolver.options);
		else
			error(['Element H{' num2str(a) ',' num2str(a) '} has absolute value larger than 1.']);
		end
	end
	E = real(Y.'*rho);
	B = reshape(B,[D,D,chi]);
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
end
end

function [B,E] = fixedblock_multicell(H,A,direction,settings)
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

if direction == 'l'
	% Calculate left block
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
	if ~all(cellfun(@(h) h{a,a} == 1,H))
		error('One or more element in H is not 1.');
	end
	B(:,a) = id;
	for a = (chi-1):-1:1
		% Compute new Y_a
		Y = applyT(reshape(B,[D,D,chi]),A,H,A,'l');
		Y = reshape(Y(:,:,a),[D^2,1]);
		% Load initial guess
		if isstruct(settings) & isfield(settings.advice,'B')
			settings.linsolver.options.v0 = reshape(settings.advice.B(:,:,a),[D*D,1]);
		end
		% Compute new L_a
		if ~all(cellfun(@(h) isempty(h{a,a}) || isscalar(h{a,a}),H))
			error('One or more diagonal elements in H are not scalar.');
		elseif any(cellfun(@(h) isempty(h{a,a}),H))
			B(:,a) = Y;
		elseif all(cellfun(@(h) h{a,a} == 1,H))
			applyM = @(x) x - applyTv(x,A,1,A,'l') + (rho.'*x)*id;
			[B(:,a),~] = linsolver(applyM,Y - (rho.'*Y)*id,settings.linsolver.options);
		elseif all(cellfun(@(h) h{a,a} < 1,H))
			applyM = @(x) x - prod(cellfun(@(h) h{a,a},H))*applyTv(x,A,1,A,'l');
			[B(:,a),~] = linsolver(applyM,Y,settings.linsolver.options);
		else
			error('One or more element in H have absolute value larger than 1.');
		end
	end
	B = reshape(B,[D,D,chi]);
	E = real(Y.'*rho)/N;
elseif direction == 'r'
	% Calculate B block
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
	if ~all(cellfun(@(h) h{a,a} == 1,H))
		error('One or more element in H is not 1.');
	end
	B(:,a) = id;
	for a = 2:chi
		% Compute new Y_a
		Y = applyT(reshape(B,[D,D,chi]),A,H,A,'r');
		Y = reshape(Y(:,:,a),[D^2,1]);
		% Load initial guess
		if isstruct(settings) & isfield(settings.advice,'B')
			settings.linsolver.options.v0 = reshape(settings.advice.B(:,:,a),[D*D,1]);
		end
		% Compute new R_a 
		if ~all(cellfun(@(h) isempty(h{a,a}) || isscalar(h{a,a}),H))
			error('One or more diagonal elements in H are not scalar.')
		elseif any(cellfun(@(h) isempty(h{a,a}),H))
			B(:,a) = Y;
		elseif all(cellfun(@(h) h{a,a} == 1,H))
			applyM = @(x) x - applyTv(x,A,1,A,'r') + (rho.'*x)*id;
			[B(:,a)] = linsolver(applyM,Y - (rho.'*Y)*id,settings.linsolver.options);
		elseif all(cellfun(@(h) h{a,a} < 1,H))
			applyM = @(x) x - H{a,a}*applyTv(x,A,1,A,'r');
			[B(:,a),~] = linsolver(applyM,Y,settings.linsolver.options);
		else
			error('One or more element in H have absolute value larger than 1.');
		end
	end
	B = reshape(B,[D,D,chi]);
	E = real(Y.'*rho)/N;
end
end
 
