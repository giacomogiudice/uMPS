import('algorithms.*');
import('util.*');
% The good ol' Pauli matrices
[X,Y,Z,I] = paulis();
% Define Ising 2-body Hamiltonian term
% Note: the constructor for TwoSiteOperator is just some eye-candy.
% Alternatively one can make the two-body gate by hand such as
% - ncon({Z,Z},{[-1,-3],[-2,-4]}) ...
% - h/2*(ncon({X,I},{[-1,-3],[-2,-4]}) + ncon({I,X},{[-1,-3],[-2,-4]}))
% Just remember to cast it as a TwoSiteOperator (otherwise it will be
% treated as an MPO).
h = 0.8;
H_twosite = TwoSiteOperator({-Z,Z},{-h,X});
% Exact energy
flambda = @(k) sqrt(h^2 + 2*h*cos(k) + 1);
E_exact = integral(@(k) (-1/pi)*flambda(k),0,pi,'RelTol',eps);
if h < 1
	mz_exact = (1 - h^2)^(1/8);
else
	mz_exact = 0;
end
mx_exact = integral(@(k) (1/pi)*((h + cos(k))./flambda(k)),0,pi,'RelTol',eps);

% Define parameters for VUMPS simulation
D = 16;
d = 2;
settings.maxit = 10;
settings.tol = eps;
% For small bond dimensions, dynamic precision can be disabled for faster convergence
% However, one should make sure that these tolerances are not smaller than settings.tol
settings.eigsolver.options.dynamictol = false;
settings.linsolver.options.dynamictol = false;
if exist('A_left','var') & exist('A_right','var') & exist('C','var')
	settings.initial.A_left = A_left;
	settings.initial.A_right = A_right;
	settings.initial.C = C;
end
% Launch VUMPS simulation
[A_left,A_right,C,~,output,B,stats] = vumps(H_twosite,D,d,settings);
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
plot(diag(C),'x')
hold on
set(gca,'yscale','log')
xlabel('$k$')
ylabel('$\lambda_k$')

% Compare observables
E = output.energy;
mx = abs(trace(applyT(C'*C,A_right,X,A_right,'l')));
mz = abs(trace(applyT(C'*C,A_right,Z,A_right,'l')));
fprintf('<H> error: %.4g, <X> error: %.4g, <Z> error: %.4g\n',abs(E - E_exact),abs(mx - mx_exact),abs(mz - mz_exact));
