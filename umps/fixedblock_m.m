function [X,E] = fixedblock_m(H,A_left,A_right,direction,settings)
N = length(H);
assert(isequal(size(H),size(A_left),size(A_right)),'Size mismatch in input cell arrays.');
chi = size(H{1},1);
[D,~,d] = size(A_left{1});
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
X = zeros([D*D,chi]);
if direction == 'l'
	% Calculate left block
	% Compute right dominant eigenvector of A_left
	% TODO add advice
	fapplyTv = @(v) applyTv(v,A_left,1,A_left,'r');
	[rho_right,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	rho_right = rho_right/(id'*rho_right);
	% Start from the corner
	a = chi;
	if ~all(cellfun(@(h) h{a,a} == 1,H))
		error('One or more element in H is not 1.');
	end
	X(:,a) = id;
	for a = (chi-1):-1:1
		% Compute new Y_a
		Y = applyT(reshape(X,[D,D,chi]),A_left,H,A_left,'l');
		Y = reshape(Y(:,:,a),[D^2,1]);
		% TODO Load initial guess
		% Compute new L_a
		if ~all(cellfun(@(h) isempty(h{a,a}) || isscalar(h{a,a}),H))
			error('One or more diagonal elements in H are not scalar.')
		elseif any(cellfun(@(h) isempty(h{a,a}),H))
			X(:,a) = Y;
		elseif all(cellfun(@(h) h{a,a} == 1,H))
			applyM_left = @(x) x - applyTv(x,A_left,1,A_left,'l') + (rho_right.'*x)*id;
			[X(:,a),~] = linsolver(applyM_left,Y - (rho_right.'*Y)*id,settings.linsolver.options);
		elseif all(cellfun(@(h) h{a,a} < 1,H))
			applyM_left = @(x) x - prod(cellfun(@(h) h{a,a},H))*applyTv(x,A_left,1,A_left,'l');
			[X(:,a),~] = linsolver(applyM_left,Y,settings.linsolver.options);
		else
			error('One or more element in H have absolute value larger than 1.')
		end
	end
	X = reshape(X,[D,D,chi]);
	E = real(Y.'*rho_right)/N;
else if direction == 'r'
	% Calculate X block
	% Compute left dominant eigenvector of A_right
	fapplyTv = @(v) applyTv(v,A_right,1,A_right,'l');
	[rho_left,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
	rho_left = rho_left/(id'*rho_left);
	% Start from the corner
	a = 1;
	if ~all(cellfun(@(h) h{a,a} == 1,H))
		error('One or more element in H is not 1.');
	end
	X(:,a) = id;
	for a = 2:chi
		% Compute new Y_a
		Y = applyT(reshape(X,[D,D,chi]),A_right,H,A_right,'r');
		Y = reshape(Y(:,:,a),[D^2,1]);
		% TODO Load initial guess
		% Compute new R_a 
		if ~all(cellfun(@(h) isempty(h{a,a}) || isscalar(h{a,a}),H))
			error('One or more diagonal elements in H are not scalar.')
		elseif any(cellfun(@(h) isempty(h{a,a}),H))
			X(:,a) = Y;
		elseif all(cellfun(@(h) h{a,a} == 1,H))
			applyM_right = @(x) x - applyTv(x,A_right,1,A_right,'r') + (rho_left.'*x)*id;
			[X(:,a)] = linsolver(applyM_right,Y - (rho_left.'*Y)*id,settings.linsolver.options);
		elseif all(cellfun(@(h) h{a,a} < 1,H))
			applyM_right = @(x) x - H{a,a}*applyTv(x,A_right,1,A_right,'r');
			[X(:,a),~] = linsolver(applyM_right,Y,settings.linsolver.options);
		else
			error(['Element H{' num2str(a) ',' num2str(a) '} has absolute value larger than 1'])
		end
	end
	X = reshape(X,[D,D,chi]);
	E = real(Y.'*rho_left)/N;
end
end