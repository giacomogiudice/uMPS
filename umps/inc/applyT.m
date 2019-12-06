function W = applyT(M,A_1,O,A_2,direction)
% APPLYT Computes the application of a tensor to the transfer matrix.
%
% W = applyT(M,A_1,O,A_2,direction)
% Computes the application of conj(A_1)^k O^kl A_2^l to M. The direction
% character specifies if on the side specified by the direction character.
% O can be a scalar, in which case the transfer matrix reduces to
% O*conj(A_1)^k A_2^k, (Be careful: if it is zero, so will be W). This
% function supports larger unit cells, by inputting A_1, A_2 as cell
% arrays.
%
% INPUT
% M         - rank-2 or rank-3 tensor for which the first index of M is
%           associated to A_1, the second to A_2, and the third to O.
% A_1,A_2   - rank-3 tensor or cell array of rank-3 tensors.
% O         - scalar, matrix, rank-3 tensor or cell array.
% direction - character specifying the side on which to apply M.
%
% OUTPUT
% W    - rank-2 or rank-3 tensor, depending on whether O is a rank-3
%        tensor.
%
% EXAMPLE
%   D = 8; d = 3; chi = 2;
%   M = randn([D,D]);
%   A = randn([D,D,d]);
%   W = applyT(M,A,1,M,'r');

% Decide which subfunction to use based on the type of the inputs
if iscell(A_1)
    W = applyT_multicell(M,A_1,O,A_2,direction);
    return
end
switch class(O)
    case 'double'
        W = applyT_generic(M,A_1,O,A_2,direction);
    case 'SchurOperator'
        W = applyT_schur(M,A_1,O,A_2,direction);
    otherwise
        error(['Unsupported function for class ' class(O)]);
end
end
