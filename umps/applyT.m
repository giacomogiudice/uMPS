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
elseif iscell(O)
    W = applyT_schur(M,A_1,O,A_2,direction);
else
    W = applyT_generic(M,A_1,O,A_2,direction);
end
end

function W = applyT_generic(M,A_1,O,A_2,direction)
if direction == 'l'
    if isempty(O)
        W = zeros([size(A_1,1),size(A_2,1)]);
    elseif isscalar(O)
        W = O*ncon({M,conj(A_1),A_2},{[1,2],[1,-1,3],[2,-2,3]});
    elseif ndims(M) == 2
        W = ncon({M,conj(A_1),O,A_2},{[1,4],[1,-1,2],[2,3],[4,-2,3]});
    else
        W = ncon({M,conj(A_1),O,A_2},{[1,5,3],[1,-1,2],[3,-3,2,4],[5,-2,4]});
    end
elseif direction == 'r'
    if isempty(O)
        W = zeros([size(A_1,2),size(A_2,2)]);
    elseif isscalar(O)
        W = O*ncon({M,conj(A_1),A_2},{[1,2],[-1,1,3],[-2,2,3]});
    elseif ndims(M) == 2
        W = ncon({M,conj(A_1),O,A_2},{[1,4],[-1,1,2],[2,3],[-2,4,3]});
    else
        W = ncon({M,conj(A_1),O,A_2},{[1,5,3],[-1,1,2],[-3,3,2,4],[-2,5,4]});
    end
else
    error(['Unrecognized direction ' direction '.']);
end
end

function W = applyT_schur(M,A_1,O,A_2,direction)
chi = size(O,1);
if direction == 'l'
    W = zeros([size(A_1,2),size(A_2,2),size(O,2)]);
    for a = chi:-1:1
        for b = chi:-1:a
            if ~isempty(O{b,a})
                    W(:,:,a) = W(:,:,a) + applyT(M(:,:,b),A_1,O{b,a},A_2,'l');
            end
        end
    end
elseif direction == 'r'
    W = zeros([size(A_1,1),size(A_2,1),size(O,1)]);
    for a = 1:chi
        for b = 1:a
            if ~isempty(O{a,b})
                    W(:,:,a) = W(:,:,a) + applyT(M(:,:,b),A_1,O{a,b},A_2,'r');
            end
        end
    end
else
    error(['Unrecognized direction ' direction '.']);
end
end

function M = applyT_multicell(M,A_1,O,A_2,direction)
N = length(A_1);
if ~iscell(O)
    O = repmat({O},1,N);
end
if direction == 'l'
    for n = 1:N
        M = applyT(M,A_1{n},O{n},A_2{n},'l');
    end
elseif direction == 'r'
    for n = N:(-1):1
        M = applyT(M,A_1{n},O{n},A_2{n},'r');
    end
else
    error(['Unrecognized direction ' direction '.']);
end
end
