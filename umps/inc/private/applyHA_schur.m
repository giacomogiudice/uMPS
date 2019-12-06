function A_new = applyHA_schur(A,H,B_left,B_right)
[D,~,d] = size(A);
chi = size(H,1);
A_new = zeros(D,D,d);
for a = 1:chi
    for b = 1:a
        if ~isempty(H{a,b})
            if isscalar(H{a,b})
                A_new = A_new + H{a,b}*ncon({B_left(:,:,a),A,B_right(:,:,b)},{[-1,1],[1,2,-3],[-2,2]});
            else
                A_new = A_new + ncon({B_left(:,:,a),A,H{a,b},B_right(:,:,b)},{[-1,1],[1,2,3],[-3,3],[-2,2]});
            end
        end
    end
end
end
