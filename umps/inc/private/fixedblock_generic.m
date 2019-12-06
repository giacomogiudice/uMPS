function [B,E] = fixedblock_generic(H,A,direction,settings)
if exist('settings','var') && isfield(settings,'advice') && isfield(settings.advice,'B')
    settings.eigsolver.options.v0 = settings.advice.B;
end
[B,E] = fixedpoint(A,H,A,direction,settings);
end
