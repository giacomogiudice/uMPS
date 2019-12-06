function [A,C_left,C_right] ...
    = solve_local(A_left,A_right,C_left,C_right,A,H,B_left,B_right,settings)
% SOLVE_LOCAL Solve local problem in VUMPS optimization.
[D,~,d] = size(A);
% Define solvers
eigsolver = settings.eigsolver.handle;
if settings.isreal && isequal(settings.eigsolver.mode,'sr')
    settings.eigsolver.mode = 'sa';
end

% Solve effective problem for A
settings.eigsolver.options.v0 = A;
fapplyHA = @(t) applyHA(t,H,B_left,B_right,A_left,A_right);
[A,~] = eigsolver(fapplyHA,[D,D,d],1,settings.eigsolver.mode,settings.eigsolver.options);
if isempty(C_right)
    % Solve effective problem for C_left
    settings.eigsolver.options.v0 = C_left;
    fapplyHC = @(t) applyHC(t,H,B_left,B_right,A_left,A_right);
    [C_left,~] = eigsolver(fapplyHC,[D,D],1,settings.eigsolver.mode,settings.eigsolver.options);
    C_left = C_left/sign(C_left(1,1));
    C_right = [];
else
    % Solve effective problem for C_left
    B_mid = applyT(B_right,A_right,H,A_right,'r');
    settings.eigsolver.options.v0 = C_left;
    fapplyHC = @(t) applyHC(t,H,B_left,B_mid);
    [C_left,~] = eigsolver(fapplyHC,[D,D],1,settings.eigsolver.mode,settings.eigsolver.options);
    C_left = C_left/sign(C_left(1,1));
    % Solve effective problem for C_right
    B_mid = applyT(B_left,A_left,H,A_left,'l');
    settings.eigsolver.options.v0 = C_right;
    fapplyHC = @(t) applyHC(t,H,B_mid,B_right);
    [C_right,~] = eigsolver(fapplyHC,[D,D],1,settings.eigsolver.mode,settings.eigsolver.options);
    C_right = C_right/sign(C_right(1,1));
end
end
