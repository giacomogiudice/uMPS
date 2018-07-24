function [B,E] = fixedblock_m(H,A,direction,settings)
N = length(H);
assert(isequal(size(H),size(A)),'Size mismatch in input cell arrays.');
chi = size(H{1},1);
[D,~,d] = size(A{1});
id = reshape(eye([D,D]),[D*D,1]);
% Shortcuts to solvers
if nargin < 4
	settings = [];
	linsolver = @(Afun,b) bicgstab(Afun,b,D*D);
	eigsolver = @eigs;
	eigsolver_options.isreal = false;
else
	tol = max(eps,settings.linsolver.options.tol);
	maxit = min(D*D,settings.linsolver.options.maxit);
	settings.linsolver.options.v0 = [];
	linsolver = settings.linsolver.handle;
	eigsolver = settings.eigsolver.handle;
	eigsolver_options.isreal = settings.isreal;
end
eigsolver_mode = 'lm';
eigsolver_options.issym = false;
B = zeros([D*D,chi]);
if direction == 'l'
	% Calculate left block
	% Compute right dominant eigenvector of A
	if isstruct(settings) & isfield(settings.advice,'C')
		C = settings.advice.C;
		eigsolver_options.v0 = reshape(C*C',[D*D,1]);
	end
	fapplyTv = @(v) applyTv(v,A,1,A,'r');
	[rho_right,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	rho_right = rho_right/(id'*rho_right);
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
			error('One or more diagonal elements in H are not scalar.')
		elseif any(cellfun(@(h) isempty(h{a,a}),H))
			B(:,a) = Y;
		elseif all(cellfun(@(h) h{a,a} == 1,H))
			applyM_left = @(x) x - applyTv(x,A,1,A,'l') + (rho_right.'*x)*id;
			[B(:,a),~] = linsolver(applyM_left,Y - (rho_right.'*Y)*id,settings.linsolver.options);
		elseif all(cellfun(@(h) h{a,a} < 1,H))
			applyM_left = @(x) x - prod(cellfun(@(h) h{a,a},H))*applyTv(x,A,1,A,'l');
			[B(:,a),~] = linsolver(applyM_left,Y,settings.linsolver.options);
		else
			error('One or more element in H have absolute value larger than 1.')
		end
	end
	B = reshape(B,[D,D,chi]);
	E = real(Y.'*rho_right)/N;
else if direction == 'r'
	% Calculate B block
	% Compute left dominant eigenvector of A
	if isstruct(settings) & isfield(settings.advice,'C')
		C = settings.advice.C;
		eigsolver_options.v0 = reshape(C*C',[D*D,1]);
	end
	fapplyTv = @(v) applyTv(v,A,1,A,'l');
	[rho_left,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	rho_left = rho_left/(id'*rho_left);
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
			applyM_right = @(x) x - applyTv(x,A,1,A,'r') + (rho_left.'*x)*id;
			[B(:,a)] = linsolver(applyM_right,Y - (rho_left.'*Y)*id,settings.linsolver.options);
		elseif all(cellfun(@(h) h{a,a} < 1,H))
			applyM_right = @(x) x - H{a,a}*applyTv(x,A,1,A,'r');
			[B(:,a),~] = linsolver(applyM_right,Y,settings.linsolver.options);
		else
			error(['Element H{' num2str(a) ',' num2str(a) '} has absolute value larger than 1'])
		end
	end
	B = reshape(B,[D,D,chi]);
	E = real(Y.'*rho_left)/N;
end
end