function [lambda] = subleading_eigenvalues(A,C,direction,n,settings)
if ~exist('n','var') || isempty(n)
	n = 1;
end
if exist('settings','var') && isfield(settings,'eigsolver')
	eigsolver = settings.eigsolver.handle;
	eigsolver_options = settings.eigsolver.options;
else
	eigsolver = linalg.eigs;
	eigsolver_options.isreal = isreal(A);
end
eigsolver_mode = 'lm';
eigsolver_options.issym = false;
eigsolver_options.fail = 'keep';

[D,~,d] = size(A);
id = eye(D,D);
if direction == 'l'
	rho = reshape((C*C').',[D*D,1]);
elseif direction == 'r'
	rho = reshape(C'*C',[D*D,1]);
else
	error(['Unrecognized direction ' direction '.']);
end

fapply = @(t) applyT(t,A,[],A,direction) - (rho(:).'*t(:))*id;
lambda = eigsolver(fapplyT,[D,D],n,eigsolver_mode,eigsolver_options);
end
