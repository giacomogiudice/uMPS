function W = applyT_generic(M,A_1,O,A_2,direction)
if direction == 'l'
    if isempty(O)
        W = ncon({M,conj(A_1),A_2},{[1,2],[1,-1,3],[2,-2,3]});
    elseif isscalar(O)
        W = O*ncon({M,conj(A_1),A_2},{[1,2],[1,-1,3],[2,-2,3]});
    elseif ndims(M) == 2
        W = ncon({M,conj(A_1),O,A_2},{[1,4],[1,-1,2],[2,3],[4,-2,3]});
    else
        W = ncon({M,conj(A_1),O,A_2},{[1,5,3],[1,-1,2],[3,-3,2,4],[5,-2,4]});
    end
elseif direction == 'r'
    if isempty(O)
        W = ncon({M,conj(A_1),A_2},{[1,2],[-1,1,3],[-2,2,3]});
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
