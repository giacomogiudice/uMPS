function nulltensor = nullspace(A,direction)
% NULLSPACE Compute the nullspace of an MPS tensor.
%
% nulltensor = nullspace(A,direction)
% Depending on the direction, computes the right or left nullspace of a
% tensor by QR decomposition.
%
% See also: QR, VUMPS.

[D,~,d] = size(A);
if direction == 'l'
    A_left = reshape(permute(A,[1,3,2]),[d*D,D]);
    [Q,~] = qr(A_left);
    nulltensor = Q(:,D+1:d*D);
    nulltensor = permute(reshape(nulltensor,[D,d,(d-1)*D]),[1,3,2]);
elseif direction == 'r'
    A_right = reshape(A,[D,d*D]);
    [Q,~] = qr(A_right.');
    nulltensor = Q(:,D+1:d*D).';
    nulltensor = reshape(nulltensor,[(d-1)*D,D,d]);
else
    error(['Unrecognized direction ' direction '.']);
end
end
