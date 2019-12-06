function [B,E] = fixedblock(H,A,direction,settings)
% FIXEDBLOCK Compute environment tensor for an MPS tensor.
%
% [B,E] = fixedblock(H,A,direction,settings) For a canonical tensor A
% (depending on the direction), computes the environment tensor B
% correspoding to the infinite contraction and the respective 'energy' E.
% The settings structure is used to to pass the initial guesses, as well
% as the internal solvers.
%
% See also: VUMPS.

if ~exist('settings','var') || isempty(settings)
    settings = struct();
end
switch class(H)
    case 'double'
        [B,E] = fixedblock_generic(H,A,direction,settings);
    case 'SchurOperator'
        [B,E] = fixedblock_schur(H,A,direction,settings);
    case 'TwoSiteOperator'
        [B,E] = fixedblock_twosite(H,A,direction,settings);
    case 'cell'
        [B,E] = fixedblock_multicell(H,A,direction,settings);
    otherwise
        error('Unrecognized class of input operator.');
end
end
