function [H_left,H_right,E] = fixedblocks_s(H,A_left,A_right,settings)
chi = size(H,1);
[D,~,d] = size(A_left);
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

% Calculate H_left block
H_left = zeros([D*D,chi]);
% Compute right dominant eigenvector of A_left
if isstruct(settings) & isfield(settings.advice,'C')
	C = settings.advice.C;
	eigsolver_options.v0 = reshape(C*C',[D*D,1]);
end
fapplyTv = @(v) applyTv(v,A_left,1,A_left,'r');
[rho_right,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
rho_right = rho_right/(id'*rho_right);
% Start from the corner
a = chi;
if H{a,a} ~= 1
	error(['Element H{' num2str(a) ',' num2str(a) '} must be 1.']);
end
H_left(:,a) = id;
for a = (chi-1):-1:1
	% Compute new Y_a
	Y = zeros([D^2,1]);
	for b = chi:-1:(a+1)
		if ~isempty(H{b,a})
			Y = Y + applyTv(H_left(:,b),A_left,H{b,a},A_left,'l');
		end
	end
	% Load initial guess
	if isstruct(settings) & isfield(settings.advice,'H_left.')
		settings.linsolver.options.v0 = reshape(settings.advice.H_left(:,:,a),[D*D,1]);
	end
	% Compute new L_a
	if ~isempty(H{a,a}) & ~isscalar(H{a,a})
		error(['Element H{' num2str(a) ',' num2str(a) '} is not a scalar.'])
	elseif isempty(H{a,a}) || H{a,a} == 0
		H_left(:,a) = Y;
	elseif H{a,a} == 1
		applyM_left = @(x) x - applyTv(x,A_left,1,A_left,'l') + (rho_right.'*x)*id;
		[H_left(:,a),~] = linsolver(applyM_left,Y - (rho_right.'*Y)*id,settings.linsolver.options);
	elseif abs(H{a,a}) < 1
		applyM_left = @(x) x - H{a,a}*applyTv(x,A_left,1,A_left,'l');
		[H_left(:,a),~] = linsolver(applyM_left,Y,settings.linsolver.options);
	else
		error(['Element H{' num2str(a) ',' num2str(a) '} has absolute value larger than 1.'])
	end
end
H_left = reshape(H_left,[D,D,chi]); 
E_left = real(Y.'*rho_right);
% Calculate H_right block
H_right = zeros([D*D,chi]);
% Compute left dominant eigenvector of A_right
if isstruct(settings) & isfield(settings.advice,'C')
	C = settings.advice.C;
	eigsolver_options.v0 = reshape(C*C',[D*D,1]);
end
fapplyTv = @(v) applyTv(v,A_right,1,A_right,'l');
[rho_left,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
rho_left = rho_left/(id'*rho_left);
% Start from the corner
a = 1;
if H{a,a} ~= 1
	error(['Element H{' num2str(a) ',' num2str(a) '} must be 1.']);
end
H_right(:,a) = id;
for a = 2:chi
	% Compute new Y_a
	Y = zeros([D^2,1]);
	for b = 1:(a-1)
		if ~isempty(H{a,b})
			Y = Y + applyTv(H_right(:,b),A_right,H{a,b},A_right,'r');
		end
	end
	% Load initial guess
	if isstruct(settings) & isfield(settings.advice,'H_right')
		settings.linsolver.options.v0 = reshape(settings.advice.H_right(:,:,a),[D*D,1]);
	end
	% Compute new R_a 
	if ~isempty(H{a,a}) & ~isscalar(H{a,a})
		error(['Element H{' num2str(a) ',' num2str(a) '} is not a scalar.'])
	elseif isempty(H{a,a}) || H{a,a} == 0
		H_right(:,a) = Y;
	elseif H{a,a} == 1
		applyM_right = @(x) x - applyTv(x,A_right,1,A_right,'r') + (rho_left.'*x)*id;
		[H_right(:,a)] = linsolver(applyM_right,Y - (rho_left.'*Y)*id,settings.linsolver.options);
	elseif abs(H{a,a}) < 1
		applyM_right = @(x) x - H{a,a}*applyTv(x,A_right,1,A_right,'r');
		[H_right(:,a),~] = linsolver(applyM_right,Y,settings.linsolver.options);
	else
		error(['Element H{' num2str(a) ',' num2str(a) '} has absolute value larger than 1.'])
	end
end
E_right = real(Y.'*rho_left);
H_right = reshape(H_right,[D,D,chi]);
E = mean([E_left,E_right]);
end