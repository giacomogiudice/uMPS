import('algorithms.*');
% The spin-1/2 operators
[sx,sy,sz] = util.su2gen(2);
% Choose model to simulate
model = 'xxz'
if isequal(model,'ising')
	% Transverse field
	h = 0.6;
	% Define MPO
	W = SchurOperator({-sz,sz},{-h,sx});
	% Exact energy
	E_exact = integral(@(k) -1/(2*pi)*sqrt(1/4+h^2+h*cos(k)),0,pi,'RelTol',eps);
elseif isequal(model,'xxz')
	% Coupling along z (set between -1 and 1)
	Delta = 2;
	% Define MPO
	W = SchurOperator({-sx,sx},{-sy,sy},{Delta,sz,sz});
	% Exact energy
	if Delta == 1
		E_exact = 1/4 - log(2);
	elseif Delta > 1
		phi = acosh(Delta);
		E_exact = Delta/4 - sinh(phi)/2 - 2*sinh(phi)*sum(1./(exp(2*phi.*(1:1e4))+1));
	else
		gamma = acos(-Delta);
		E_exact = Delta/4 - sin(gamma)*integral(@(x)(1-tanh(x*gamma)./tanh(x*pi)),0,inf,'RelTol',eps);
	end
end
% A 4-site unit cell is overkill, only a 2-site is necessary
H = {W,W,W,W};

% Define parameters for VUMPS simulation
D = [8,16,24,32];
d = 2;
settings.maxit = 20;
settings.tol = 1e-12;
if exist('A_left','var') & exist('A_right','var') & exist('C','var')
	settings.initial.A_left = A_left;
	settings.initial.A_right = A_right;
	settings.initial.C = C;
end
% Launch VUMPS simulation
[A_left,A_right,C,~,output,~,stats] = vumps(H,D,d,settings);
output

% Plot results
figure(1)
plot(1:output.iter,abs(stats.energy - E_exact),'-o')
hold on
plot(1:output.iter,abs(stats.energydiff),'-+')
plot(1:output.iter,stats.err,'-s')
set(gca,'yscale','log')
xlabel('iterations')
legend({'$|E - E_{\rm exact}|$','$|E^{(n)} - E^{(n-1)}|$','$\epsilon$'})

figure(2)
for n = 1:length(H)
	plot(diag(C{n}),'x')
	hold on
end
set(gca,'yscale','log')
xlabel('$k$')
ylabel('$\lambda_k$')
