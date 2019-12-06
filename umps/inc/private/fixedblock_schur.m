function [B,E] = fixedblock_schur(H,A,direction,settings)
import('linalg.geometric_series');
chi = size(H,1);
D = size(A,1);
B = zeros([D,D,chi]);

% Compute dominant eigenvectors of A
[rho_left,rho_right] = fixedpoint_canonical(A,direction,settings);

% Start from corner and solver for each B_a
if direction == 'l'
    corner = chi;
    ind = (chi-1):-1:1;
    rho = rho_right;
else
    corner = 1;
    ind = 2:1:chi;
    rho = rho_left;
end
assert(H{corner,corner} == 1,'Element H{%d,%d} must be 1.',corner,corner);
B(:,:,corner) = eye(D,D);
for a = ind
    % Compute new Y_a
    Y = zeros([D,D]);
    if direction == 'l'
        for b = chi:-1:(a+1)
            if ~isempty(H{b,a})
                Y = Y + applyT(B(:,:,b),A,H{b,a},A,'l');
            end
        end
    else
        for b = 1:(a-1)
            if ~isempty(H{a,b})
                Y = Y + applyT(B(:,:,b),A,H{a,b},A,'r');
            end
        end
    end
    % Compute new B_a
    if isempty(H{a,a})
        B(:,:,a) = Y;
    else
        assert(isscalar(H{a,a}),'Element H{%d,%d} is not a scalar.',a,a);
        % Load initial guess
        if isstruct(settings) && isfield(settings.advice,'B')
            settings.linsolver.options.v0 = settings.advice.B(:,:,a);
        end
        if H{a,a} == 1
            v_left = rho_left;
            v_right = rho_right;
        elseif abs(H{a,a}) < 1
            v_left = 0;
            v_right = 0;
        else abs(H{a,a}) > 1
            error('Diagonal element in H has absolute value larger than 1.');
        end
        Tfun = @(x) applyT(x,A,H{a,a},A,direction);
        B(:,:,a) = geometric_series(Y,Tfun,v_left,v_right,direction,settings);
    end
end
E = real(Y(:).'*rho(:));
end
