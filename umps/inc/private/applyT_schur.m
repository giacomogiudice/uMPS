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
