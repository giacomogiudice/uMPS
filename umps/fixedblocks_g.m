function [L,R,E] = fixedblocks_g(H,A_left,A_right,settings)
chi = size(H,1);
[D,~,d] = size(A_left);
if nargin < 4
	settings = [];
	eigsolver = @eigs;
	eigsolver_options.isreal = false;
else
	eigsolver = settings.eigsolver.handle;
	eigsolver_options.isreal = settings.isreal;
end
eigsolver_mode = 'lm';
eigsolver_options.issym = true;

% Compute right dominant eigenvector
fapplyTv = @(v) applyTv(v,A_left,H,A_left,'l');
[L,eta_left] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
L = reshape(L,[D,D,chi]);

% Compute right dominant eigenvector
fapplyTv = @(v) applyTv(v,A_right,H,A_right,'r');
[R,eta_right] = eigsolver(fapplyTv,D^2,1,eigsolver_mode,eigsolver_options);
R = reshape(R,[D,D,chi]);
E = mean([eta_left eta_right]);
end