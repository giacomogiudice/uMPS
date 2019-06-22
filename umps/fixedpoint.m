function [rho,eta,rho_adj] = fixedpoint(A,direction)
% Finds right fixedpoint by default
rho_adj = [];
if nargin < 2
	direction = 'r';
elseif strcmp(direction,'all')
	[rho,eta] = fixedpoint(A,'l');
	[rho_adj,eta_adj] = fixedpoint(A,'r');
	rho_norm = sqrt(trace(rho*rho_adj.'));
	rho = rho/rho_norm;
	rho_adj = rho_adj/rho_norm;
	eta = mean([eta,eta_adj]);
	return
end
D = size(A,1);
options.isreal = isreal(A);
options.issym = 0;
[rho,eta] = eigs(@(v) applyTv(v,A,1,A,direction),D^2,1,'lm',options);
rho = reshape(rho,[D,D]);
rho = rho/trace(rho);
end
