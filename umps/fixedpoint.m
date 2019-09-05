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
% See also: APPLYT, VUMPS.

if exist('settings','var') && isfield(settings,'eigsolver')
    eigsolver = settings.eigsolver.handle;
    eigsolver_options = settings.eigsolver.options;
else
    eigsolver = @eigs;
    eigsolver_options.isreal = isequal(A_1,A_2) && all([isreal(A_1),isreal(O)]);
end
eigsolver_mode = 'lm';
eigsolver_options.issym = false;
eigsolver_options.fail = 'keep';
if iscell(A_1) && iscell(A_2)
    D_1 = size(A_1{1},1);
    D_2 = size(A_2{1},1);
else
    D_1 = size(A_1,1);
    D_2 = size(A_2,1);
end
if ndims(O) == 4
    chi = size(O,1);
else
    chi = 1;
end
fun = @(v) applyTv(v,A_1,O,A_2,direction);
if iscell(A_1) && iscell(O) && ~iscell(O{1})
    fun = @(v) applyTvNtimes(v,A_1,H,A_2,direction);
end
[rho,eta] = eigsolver(fun,D_1*D_2*chi,1,eigsolver_mode,eigsolver_options);
rho = reshape(rho,[D_1,D_2,chi]);
end

function v = applyTvNtimes(v,A1,H,A2,direction)
N = length(A1);
if direction == 'l'
    ind = 1:N;
elseif direction == 'r'
    ind = N:(-1):1;
end
for n = ind
    v = applyTv(v,A1{n},H{n},A2{n},direction);
end
end
