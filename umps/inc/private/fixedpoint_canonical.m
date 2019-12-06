function [rho_left,rho_right] = get_fixedpoint_canonical(A,direction,settings)
if iscell(A)
    D = size(A{1},1);
else
    D = size(A,1);
end
id = eye(D,D);
if direction == 'l'
    % Left block
    if isstruct(settings) && isfield(settings.advice,'C')
        C = settings.advice.C;
        settings.eigsolver.options.v0 = (C*C').';
    end
    rho = fixedpoint(A,[],A,'r',settings);
    rho_left = id;
    rho_right = rho/trace(rho);
elseif direction == 'r'
    % Right block
    if isstruct(settings) && isfield(settings.advice,'C')
        C = settings.advice.C;
        settings.eigsolver.options.v0 = C'*C;
    end
    rho = fixedpoint(A,[],A,'l',settings);
    rho_left = rho/trace(rho);
    rho_right = id;
else
    error(['Unrecognized direction ' direction '.']);
end
end
