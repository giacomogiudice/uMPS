clear;
import('algorithms.*');
import('util.*');

% Parameters for quench
h0 = 1.5;
h = 2;

% Get Hamiltonian for h0
[X,Y,Z,I] = paulis();
H = SchurOperator({-X,X},{-h0,Z});

% Compute ground state at h0
D = 32;
d = 2;
fprintf('Optimizing tensors for h = %g...\n',h0);
[A0_left,A0_right,C0,A0] = vumps(H,D,d,1e-14);

% Get Hamiltonian for h
H = SchurOperator({-X,X},{-h,Z});

% Launch time evolution
settings.niter = 40;
settings.timestep = -1i*0.1;
settings.customfun = @(A_left,A_right,C,varargin) real(trace(applyT(C'*C,A_right,Z,A_right,'l')));
fprintf('Running time-evolution for h = %g...\n',h);
[A_left,A_right,C,A,output,~,stats] = tdvp(H,A0_left,A0_right,C0,A0,settings);
mags = cell2mat(stats.customvalue);

% Compare with exact results (formula from G. Giudici's master thesis)
dt = abs(settings.timestep);
t_final = abs(dt*settings.niter);
t_sampling = linspace(0,t_final,200);
mags_exact = zeros(size(mags));
disp_rel = @(p,g) 2*sqrt(1 + g^2 - 2.*g.*cos(p));
integrand = @(p,t) ...
    ((h*h0 + 1 - (h + h0).*cos(p)).*(h - cos(p)) - (h - h0).*sin(p).^2.*cos(2*disp_rel(p,h)*t)) ...
    ./(disp_rel(p,h ).^2.*disp_rel(p,h0));
for ind = 1:length(t_sampling)
	mags_exact(ind) = 4/pi*integral(@(p) integrand(p,t_sampling(ind)),0,2*pi);
end

figure(1)
plot(dt*(1:settings.niter),mags,'x');
hold on
plot(t_sampling,mags_exact,'-');
xlabel('$t$')
ylabel('$\langle Z \rangle$')
legend('TDVP','Exact')
