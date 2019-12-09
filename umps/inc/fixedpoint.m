function [rho,eta] = fixedpoint(A_1,O,A_2,direction,settings)
% FIXEDPOINT Computes the leading eigenvector of a transfer matrix.
%
% rho = fixedpoint(A_1,O,A_2,direction)
% Computes the left ('l') or right ('r') eigenvector of the fixedpoint rho
% depending on the direction flag. The transfer matrix given by
% conj(A_1)^k O^kl A_2^l is not constructed explicitly, but just it's
% application to a vector is computed using applyT.
%
% rho = fixedpoint(A_1,O,A_2,direction,settings)
% Additionally, one can provide a VUMPS settings structure to use a
% specific eigensolver.
%
% [rho,eta] = fixedpoint(A_1,O,A_2,direction,...)
% The second output eta is the corresponding eigenvalue.
%
% See also: APPLYT.
if exist('settings','var') && isfield(settings,'eigsolver')
    eigsolver = settings.eigsolver.handle;
    eigsolver_options = settings.eigsolver.options;
else
    eigsolver = @linalg.eigs;
    eigsolver_options.isreal = isequal(A_1,A_2) && all([isreal(A_1),isreal(O)]);
end
eigsolver_mode = 'lm';
eigsolver_options.issym = false;
eigsolver_options.fail = 'keep';

if iscell(A_1) && iscell(A_2)
    % Multicell case
    D_1 = size(A_1{1},1);
    D_2 = size(A_2{1},1);
else
    % Regular case
    D_1 = size(A_1,1);
    D_2 = size(A_2,1);
end
if ndims(O) == 4
    chi = size(O,1);
else
    chi = 1;
end
fapply = @(t) applyT(t,A_1,O,A_2,direction);
[rho,eta] = eigsolver(fapply,[D_1,D_2,chi],1,eigsolver_mode,eigsolver_options);
end
