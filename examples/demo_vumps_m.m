% The good ol' Pauli matrices
sx = [0,1;1,0];
sy = [0,-1i;1i,0];
sz = [1,0;0,-1];
si = eye(2);
% Define Ising 2-body Hamiltonian term
h = 0.2;
W = cell(3,3);
W{1,1} = 1;
W{2,1} = sz;
W{3,1} = -h*sx;
W{3,2} = sz;
W{3,3} = 1;
H = {W,W};

% Exact energy
flambda = @(k) sqrt(h^2 + 2*h*cos(k) + 1);
E_exact = integral(@(k) (-1/pi)*flambda(k),0,pi,'RelTol',eps);

% Define parameters for VUMPS simulation
D = 10;
d = 2;
settings.mode = 'multicell';
settings.maxit = 20;
settings.tol = eps;

% Launch VUMPS simulation
[A_left,A_right,C,output,stats] = vumps(H,D,d,settings);
output
E = output.energy;

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
