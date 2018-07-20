function [H_left,H_right,E] = fixedblocks_twosite(h,A_left,A_right,settings)
[D,~,d] = size(A_left);
id = reshape(eye(D,D),[D*D,1]);
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

% Left block
% Compute right dominant eigenvector of A_left
if isstruct(settings) & isfield(settings.advice,'C')
	C = settings.advice.C;
	eigsolver_options.v0 = reshape(C*C',[D*D,1]);
end
fapplyTv = @(v) applyTv(v,A_left,1,A_left,'r');
[rho_right,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
rho_right = rho_right/(id'*rho_right);
% Compute left boundary
block_left = ncon({A_left,A_left},{[-1,1,-2],[1,-4,-3]});
h_left = ncon({conj(block_left),h,block_left},{[5,1,2,-1],[1,2,3,4],[5,3,4,-2]});
h_left = reshape((h_left+h_left')/2,[D*D,1]);
% Solve left linear system
if isstruct(settings) & isfield(settings.advice,'H_left')
	settings.linsolver.options.v0 = reshape(settings.advice.H_left,[D*D,1]);
end
applyM_left = @(x) x - applyTv(x,A_left,1,A_left,'l') + (rho_right.'*x)*id;
B_left = h_left - real(rho_right.'*h_left)*id;
[H_left,~] = linsolver(applyM_left,B_left,settings.linsolver.options);
H_left = reshape(H_left,[D,D]);
E_left = real(h_left.'*rho_right);

% Right block
% Compute left dominant eigenvector of A_right
if isstruct(settings) & isfield(settings.advice,'C')
	C = settings.advice.C;
	eigsolver_options.v0 = reshape(C*C',[D*D,1]);
end
fapplyTv = @(v) applyTv(v,A_right,1,A_right,'l');
[rho_left,~] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
rho_left = rho_left/(id'*rho_left);
% Compute right boundary
block_right = ncon({A_right,A_right},{[-1,1,-2],[1,-4,-3]});
h_right = ncon({conj(block_right),h,block_right},{[-1,1,2,5],[1,2,3,4],[-2,3,4,5]});
h_right = reshape((h_right+h_right)/2,[D*D,1]);
% Solve right linear system
if isstruct(settings) & isfield(settings.advice,'H_right')
	settings.linsolver.options.v0 = reshape(settings.advice.H_right,[D*D,1]);
end
applyM_right = @(x) x - applyTv(x,A_right,1,A_right,'r') + (rho_left.'*x)*id;
B_right = h_right - real(rho_left.'*h_right)*id;
[H_right,~] = linsolver(applyM_right,B_right,settings.linsolver.options);
H_right = reshape(H_right,[D,D]);
E_right = real(rho_left.'*h_right);
% Compute energy
E = mean([E_left,E_right]);
end
 