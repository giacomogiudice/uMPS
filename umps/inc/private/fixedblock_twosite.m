function [B,E] = fixedblock_twosite(H,A,direction,settings)
import('linalg.geometric_series');
D = size(A,1);
% Compute dominant eigenvectors of A
[rho_left,rho_right] = fixedpoint_canonical(A,direction,settings);

% Compute relevant boundary
block = ncon({A,A},{[-1,1,-2],[1,-4,-3]});
if direction == 'l'
    h = ncon({conj(block),H,block},{[5,1,2,-1],[1,2,3,4],[5,3,4,-2]});
    rho = rho_right;
elseif direction == 'r'
    h = ncon({conj(block),H,block},{[-1,1,2,5],[1,2,3,4],[-2,3,4,5]});
    rho = rho_left;
else
    error(['Unrecognized direction ' direction '.']);
end

% Solve linear system
if isstruct(settings) && isfield(settings.advice,'B')
    settings.linsolver.options.v0 = settings.advice.B;
end
Tfun = @(x) applyT(x,A,1,A,direction);
B = geometric_series(h,Tfun,rho_left,rho_right,direction,settings);
E = real(h(:).'*rho(:));
end
