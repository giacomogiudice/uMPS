% The spin-1/2 operators
sx = 1/2*[0,1;1,0];
sy = 1/2*[0,-1i;1i,0];
sz = 1/2*[1,0;0,-1];
% Choose model to simulate
model = 'xxz'
if isequal(model,'ising')
	% Transverse field
	h = 0.6;
	% Define MPO
	W = cell(3,3);
	W{1,1} = 1;
	W{2,1} = sz;
	W{3,1} = h*sx;
	W{3,2} = sz;
	W{3,3} = 1;
	% Exact energy
	E_exact = integral(@(k) -1/(2*pi)*sqrt(1/4+h^2+h*cos(k)),0,pi,'RelTol',eps);
elseif isequal(model,'xxz')
	% Coupling along z (set between -1 and 1)
	Delta = 2;
	% Define MPO
	W = cell(5,5);
	W{1,1} = 1;
	W{2,1} = sz;
	W{3,1} = sy;
	W{4,1} = sx;
	W{5,2} = Delta*sz;
	W{5,3} = sy;
	W{5,4} = sx;
	W{5,5} = 1;
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
H = {W,W,W,W};

% Define parameters for VUMPS simulation
D = 9;
d = 2;
settings.mode = 'multicell';
settings.maxit = 20;
settings.tol = 1e-12;
if exist('A_left','var') & exist('A_right','var') & exist('C','var')
	settings.initial.A_left = A_left;
	settings.initial.A_right = A_right;
	settings.initial.C = C;
end
% Launch VUMPS simulation
[~,~,C,~,output,~,stats] = vumps(H,D,d,settings);
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
