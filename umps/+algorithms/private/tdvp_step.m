function [A,C] = tdvp_step(dt,A_left,A_right,C_left,C_right,A,H,B_left,B_right,settings)
[D,~,d] = size(A);
% Define integrator
integrator = settings.integrator.handle;
integrator_options = settings.integrator.options;

% Perform integration step for A
fapplyHA = @(x) applyHA(x,H,B_left,B_right,A_left,A_right);
A = integrator(dt,fapplyHA,A,integrator_options);
A = A/sqrt(A(:)'*A(:));

% Solve integration step for C
fapplyHC = @(x) applyHC(x,H,B_left,B_right,A_left,A_right);
C = integrator(dt,fapplyHC,C_left,integrator_options);
C = C/sqrt(C(:)'*C(:));
end
