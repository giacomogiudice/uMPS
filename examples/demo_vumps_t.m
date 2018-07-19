% The good ol' Pauli matrices
sx = [0,1;1,0];
sy = [0,-1i;1i,0];
sz = [1,0;0,-1];
si = eye(2);
% Define Ising 2-body Hamiltonian term
g = 0.9;
h = -ncon({sz,sz},{[-1,-3],[-2,-4]}) - g/2*(ncon({sx,si},{[-1,-3],[-2,-4]}) + ncon({si,sx},{[-1,-3],[-2,-4]}));
% Exact energy
E_exact = integral(@(k) (-1/pi)*sqrt(1+g^2-2*g*cos(k)),0,pi,'RelTol',eps);

% Define parameters for VUMPS simulation
D = 20;
d = 2;
settings.mode = 'twosite';
settings.maxit = 15;
settings.tol = eps;
if exist('A_left','var') & exist('A_right','var') & exist('C','var')
	settings.initial.A_left = A_left;
	settings.initial.A_right = A_right;
	settings.initial.C = C;
end
% Launch VUMPS simulation
[A_left,A_right,C,output,stats] = vumps(h,D,d,settings);
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
set(gca,'yscale','log')
xlabel('$k$')
ylabel('$\lambda_k$')

% Compute observables
magn_x_twosite = 1/2*(ncon({sx,si},{[-1,-3],[-2,-4]}) + ncon({si,sx},{[-1,-3],[-2,-4]}));
magn_z_twosite = 1/2*(ncon({sz,si},{[-1,-3],[-2,-4]}) + ncon({si,sz},{[-1,-3],[-2,-4]}));
block_twosite = ncon({A_left,C,A_right},{[-1,1,-2],[1,2],[2,-4,-3]});
E = real(ncon({conj(block_twosite),h,block_twosite},{[5,1,2,6],[1,2,3,4],[5,3,4,6]}));
magn_x = real(ncon({conj(block_twosite),magn_x_twosite,block_twosite},{[5,1,2,6],[1,2,3,4],[5,3,4,6]}));
magn_z = real(ncon({conj(block_twosite),magn_z_twosite,block_twosite},{[5,1,2,6],[1,2,3,4],[5,3,4,6]}));
[magn_x magn_z]
