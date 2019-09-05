function y = geometric_series(x,T,v_left,v_right,direction,settings)
% GEOMETRIC_SERIES Compute the geometric sum of a matrix applied to a vector.
%
% y = geometric_series(x,T,v_left,v_right,direction)
% Depending on the direction ('l' or 'r') specified, this corresponds to
%   |y> = (1 + T + T^2 + ...)|x>    (right)
%   <y| = <x|(1 + T + T^2 + ...)    (left)
% with a divergent contribution, i.e. a unique leading eigenvalue
% <v_left|T|v_right> = 1. This is achieved by solving the linear systems
%   (1 - T + |v_left><v_right|)|y> = (1 - |v_left><v_right|)|x>    (left)
%   <y|(1 - T + |v_left><v_right|) = <x|(1 - |v_left><v_right|)    (right)
%
% y = geometric_series(x,T,v_left,v_right,direction,settings)
% Does the same as above, except using the linear solver specified in the
% settings field.
%
% INPUT
% x                 - vector to apply to the geometric series.
% T                 - matrix or function handle corresponding to the matrix
%                   in the summation.
% v_left,v_right    - dominant left and right eigenvectors of the matrix
%                   (can be empty or 0).
% direction         - specifies side on which to apply 
% settings          - (optional) settings structure for linear solver.
%
% OUTPUT
% y     - resulting output vector.

if isnumeric(T)
    if direction == 'l'
        Tfun = @(v) T.'*v;
    elseif direction == 'r'
         Tfun = @(v) T*v;
    else
    error(['Unrecognized direction ' direction '.']);
    end
else
    Tfun = T;
end

% Define linear solver
if exist('settings','var') && isfield(settings,'linsolver')
    linsolver = settings.linsolver.handle;
    linsolver_options = settings.linsolver.options;
else
    linsolver = @bicgstab_;
    linsolver_options = struct();
end

if isempty(v_left) || isscalar(v_left)
    v_left = 0;
end
if isempty(v_right) || isscalar(v_right)
    v_right = 0;
end

if direction == 'l'
    applyM = @(v) v - Tfun(v) + (v_right.'*v)*v_left;
    b = x - (v_right.'*x)*v_left;
elseif direction == 'r'
    applyM = @(v) v - Tfun(v) + (v_left.'*v)*v_right;
    b = x - (v_left.'*x)*v_right;
else
    error(['Unrecognized direction ' direction '.']);
end
[y,~] = linsolver(applyM,b,linsolver_options);
end
