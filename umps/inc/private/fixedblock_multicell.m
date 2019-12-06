function [B,E] = fixedblock_multicell(H,A,direction,settings)
import('linalg.geometric_series');
if isnumeric(H{1})
    N = length(H);
    [B,E] = fixedblock_generic(H,A,direction,settings);
    E = E/N;
    return
end
N = length(H);
assert(isequal(size(H),size(A)),'Size mismatch in input cell arrays.');
chi = size(H{1},1);
D = size(A{1},1);
B = zeros([D,D,chi]);
Hall = cell2schur(H);

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
assert(all(cellfun(@(h) h == 1,Hall{corner,corner})),'All {%d,%d} H must be 1.',corner,corner);
B(:,:,corner) = eye(D,D);
for a = ind
    % Compute new Y_a
    Y = zeros([D,D]);
    if direction == 'l'
        for b = chi:-1:a
            if ~isempty(Hall{b,a})
                for s = 1:size(Hall{b,a},1)
                    Y = Y + applyT(B(:,:,b),A,Hall{b,a}(s,:),A,'l');
                end
            end
        end
    else
        for b = 1:(a-1)
            if ~isempty(Hall{a,b})
                for s = 1:size(Hall{a,b},1)
                    Y = Y + applyT(B(:,:,b),A,Hall{a,b}(s,:),A,'r');
                end
            end
        end
    end
    % Compute new B_a
    Hdiag = Hall{a,a};
    if isempty(Hdiag)
        B(:,:,a) = Y;
    else
        assert(all(cellfun(@(h) isscalar(h),Hdiag)),'One or more diagonal elements in H are not scalar.');
        % Load initial guess
        if isstruct(settings) && isfield(settings.advice,'B')
            settings.linsolver.options.v0 = settings.advice.B(:,:,a);
        end
        if all(cellfun(@(h) h == 1,Hdiag))
            v_left = rho_left;
            v_right = rho_right;
        elseif abs(H{a,a}) < 1
            v_left = 0;
            v_right = 0;
        else
            error('Diagonal element in H has absolute value larger than 1.');
        end
        Tfun = @(x) prod(cell2mat(Hdiag))*applyT(x,A,[],A,direction);
        B(:,:,a) = geometric_series(Y,Tfun,v_left,v_right,direction,settings);
    end
end
E = real(Y(:).'*rho(:))/N;
end
