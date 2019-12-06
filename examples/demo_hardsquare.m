% Compute the hard-square constant, by computing the classical partition
% function of the hard square model.
import('algorithms.*');
% The hard square MPO (non-symmetrized)
d = 2;
G = [1,1;1,0];
X = zeros([d,d,d,d]);
X(1,1,1,1) = 1;
X(2,2,2,2) = 1;
T = ncon({G,X,G},{[1,-2],[1,-1,2,-3],[2,-4]});

% Numerical solution of the partition function by Baxter
E_exact = 1.503048082475332264322066329475553689385781;

% Define parameters for VUMPS simulation
D = 8;
settings.isreal = 1;
settings.eigsolver.mode = 'lm';
settings.maxit = 20;
settings.tol = eps;

% Launch VUMPS simulation
[A_left,A_right,C,A,output,~,stats] = vumps(T,D,d,settings);
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

