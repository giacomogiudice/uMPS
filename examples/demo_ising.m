% Solve the 2D classical Ising model by computing the fixed point of the
% partition function as an MPS.
import('algorithms.*');
% Inverse temperature parameter
beta = 0.45;

% Define the MPO associated to the transfer matrix (symmetrized)
d = 2;
Q = exp(-beta*[-1,1;1,-1])^(1/2);
X = zeros([d,d,d,d]);
X(1,1,1,1) = 1;
X(2,2,2,2) = 1;
T = ncon({Q,Q,X,Q,Q},{[-1,1],[2,-2],[1,2,3,4],[-3,3],[4,-4]});

% Here are the exact values (following Baxter's notation)
beta_crit = log(1 + sqrt(2))/2;
k_baxter = 1/(sinh(2*beta))^2;
F = @(theta) log(2*(cosh(2*beta)^2 + (1/k_baxter)*sqrt(1 + k_baxter^2 - 2*k_baxter*cos(2*theta))));
f_exact = (-1/beta)*integral(F,0,pi)/(2*pi);
m_exact = (beta > beta_crit)*(1 - k_baxter^2)^(1/8);

% Define parameters for VUMPS simulation
D = [4,8,16,32,64];	% Variable bond dimension
settings.eigsolver.mode = 'lm';
settings.maxit = 20;
settings.tol = 1e-14;
settings.isreal = true;

% Launch VUMPS simulation
[A_left,A_right,C,A,output,B] = vumps(T,D,d,settings);

% Free enery per site
f = -log(output.energy)/beta;

sz = [1,0;0,-1];
Y = ncon({X,sz},{[1,-2,-3,-4],[-1,1]});
M = ncon({Q,Q,Y,Q,Q},{[-1,1],[2,-2],[1,2,3,4],[-3,3],[4,-4]});

magn_onesite = ncon({applyHA(A,M,B.left,B.right),A},{[1,2,3],[1,2,3]});
norm_onesite = ncon({applyHA(A,T,B.left,B.right),A},{[1,2,3],[1,2,3]});
m = abs(magn_onesite/norm_onesite);

fprintf('Results for beta = %.4g, %.4g from critical point\n',beta,beta - beta_crit);
fprintf('Free energy (error): %.4g (%.4g), magnetization (error): %.4g (%.4g)\n',f,abs(f - f_exact),m,abs(m - m_exact));
